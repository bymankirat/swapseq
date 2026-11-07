#!/usr/bin/env node
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mankirat Singh
// SwapSeq on the command line: the same engine as the web page and the R package, for scripts and pipelines.
// Node.js 18 or newer, no other software. Run `swapseq --help`.
import fs from 'fs';import path from 'path';
import {analyse,annotate,distances,markerDistances,measurePower,findShifts,hostCheck,addHost} from '../method/core.mjs';

const VERSION='1.2.0';
const HELP=`SwapSeq ${VERSION}: finds swapped and mislabelled samples in longitudinal microbiome studies.

Usage:
  swapseq <feature-table> <metadata> --subject <column> [options]

Inputs:
  feature-table   Features in rows, samples in columns: TSV or CSV (QIIME 2 / biom exports work),
                  or a BIOM 1.0 JSON file. For BIOM 2 (HDF5) or .qza files, export first:
                  qiime tools export --input-path table.qza --output-path out
                  biom convert -i out/feature-table.biom -o table.tsv --to-tsv
  metadata        TSV or CSV, sample IDs in the first column.

Options:
  --subject COL       Metadata column that says who each sample came from (required)
  --type TYPE         counts (default: taxa, ASV or OTU counts) or markers (strain markers, genes or SNVs;
                      any value above 0 = present)
  --host FILE         Host-DNA read counts at common human SNPs: TSV with columns sample, site, alt, total.
                      Matches samples to people by their own DNA (newborns, antibiotics). Human genetic data:
                      keep it local and follow your ethics approval.
  --well COL, --plate COL, --household COL
                      Optional metadata columns: plate positions (finds sample-sheet shifts) and households
                      (people who live together are not over-flagged)
  --prior P           Expected share of mislabelled samples (default 0.01)
  --min-depth N       Leave out samples with fewer than N reads (counts mode; default 0 = keep all).
                      Flags on samples under 1,000 reads are marked as low depth either way.
  --samples-in-rows   The table has samples in rows and features in columns
  --out PREFIX        Write PREFIX.tsv (one row per sample) and PREFIX.json (summary). Default: swapseq
  --no-power          Skip the measured-power step (faster on very large tables)
  --fail-on-flags     Exit with code 2 when any sample is strongly flagged (for pipelines)
  --quiet             Print nothing but errors
  --version, --help

A flag is a reason to check lab records, not proof. Research use only; not for diagnosis.
Cite: Singh M. (2026) SwapSeq. https://doi.org/10.5281/zenodo.23963612`;

// ---------- arguments ----------
const argv=process.argv.slice(2),opt={type:'counts',prior:.01,out:'swapseq'},pos=[];
for(let i=0;i<argv.length;i++){const a=argv[i];
  if(a==='--help'||a==='-h'){console.log(HELP);process.exit(0)}
  if(a==='--version'||a==='-v'){console.log(VERSION);process.exit(0)}
  if(a==='--samples-in-rows'){opt.rowsAreSamples=true;continue}
  if(a==='--no-power'){opt.noPower=true;continue}
  if(a==='--fail-on-flags'){opt.failOnFlags=true;continue}
  if(a==='--quiet'){opt.quiet=true;continue}
  if(a.startsWith('--')){const k=a.slice(2),v=argv[++i];if(v==null)die(`${a} needs a value`);
    if(!['subject','type','host','well','plate','household','prior','out','min-depth'].includes(k))die(`unknown option ${a}`);opt[k]=k==='prior'||k==='min-depth'?+v:v;continue}
  pos.push(a)}
function die(m){console.error('swapseq: '+m+'\nRun swapseq --help for usage.');process.exit(1)}
const say=(...m)=>{if(!opt.quiet)console.log(...m)};
if(pos.length!==2)die('give a feature table and a metadata file');
if(!opt.subject)die('--subject is required');
if(!['counts','markers'].includes(opt.type))die('--type must be counts or markers');
if(!(opt.prior>0&&opt.prior<.5))die('--prior must be between 0 and 0.5');

// ---------- reading files ----------
const read=f=>{if(!fs.existsSync(f))die(`file not found: ${f}`);return fs.readFileSync(f,'utf8').replace(/^﻿/,'')};
const sepOf=(f,line)=>/\.csv$/i.test(f)?',':(line.split('\t').length>=line.split(',').length?'\t':',');
function splitLine(l,sep){if(sep==='\t')return l.split('\t').map(s=>s.trim());const o=[];let cur='',q=false;for(const ch of l){if(ch==='"'){q=!q;continue}if(ch===sep&&!q){o.push(cur.trim());cur=''}else cur+=ch}o.push(cur.trim());return o}
function readTable(f){
  const txt=read(f);
  if(txt.trimStart().startsWith('{')){const j=JSON.parse(txt);if(!j.rows||!j.columns||!j.data)die('this JSON file is not a BIOM 1.0 table');
    const feats=j.rows.map(r=>r.id),per=j.columns.map(()=>[]);
    if(j.matrix_type==='dense')j.data.forEach((row,r)=>row.forEach((v,c)=>{if(v>0)per[c].push([r,v])}));else for(const [r,c,v] of j.data)if(v>0)per[c].push([r,v]);
    return {feats,samples:j.columns.map((col,c)=>{const a=per[c].sort((x,y)=>x[0]-y[0]);return {id:String(col.id),idx:Int32Array.from(a.map(x=>x[0])),val:Float64Array.from(a.map(x=>x[1]))}})}}
  let L=txt.split(/\r?\n/).filter(l=>l.trim()!=='');
  L=L.filter((l,i)=>!(l.startsWith('# ')&&!/^#\s*(OTU|ASV|Feature|feature)[ _]?ID/i.test(l)));   // biom "# Constructed from" lines
  const sep=sepOf(f,L[0]),rows=L.map(l=>splitLine(l,sep)),head=rows[0];
  let feats,ids,get;
  if(!opt.rowsAreSamples){feats=rows.slice(1).map(r=>r[0]);ids=head.slice(1);get=(s,fi)=>rows[fi+1][s+1]}
  else{ids=rows.slice(1).map(r=>r[0]);feats=head.slice(1);get=(s,fi)=>rows[s+1][fi+1]}
  let bad=0;const samples=ids.map((id,s)=>{const I=[],V=[];for(let fi=0;fi<feats.length;fi++){const raw=get(s,fi);if(raw==null||raw==='')continue;const v=+raw;if(!Number.isFinite(v)){bad++;continue}if(v<0)die('the table has negative values');if(v>0){I.push(fi);V.push(v)}}
    return {id:String(id),idx:Int32Array.from(I),val:Float64Array.from(V)}});
  if(bad>feats.length)die('the table has non-numeric values; is it features in rows and samples in columns? (or add --samples-in-rows)');
  return {feats,samples}}
function readMeta(f){const L=read(f).split(/\r?\n/).filter(l=>l.trim()!=='');const sep=sepOf(f,L[0]);
  const head=splitLine(L[0],sep),out=new Map();
  for(const l of L.slice(1)){if(l.startsWith('#'))continue;const c=splitLine(l,sep),o={};head.forEach((h,i)=>o[h]=c[i]??'');out.set(c[0],o)}
  return {head,rows:out}}
function readHost(f,ids){const L=read(f).split(/\r?\n/).filter(l=>l.trim()!=='');const sep=sepOf(f,L[0]),h=splitLine(L[0],sep).map(x=>x.toLowerCase());
  const ix=['sample','site','alt','total'].map(k=>h.indexOf(k));if(ix.some(i=>i<0))die('--host file needs the columns sample, site, alt and total');
  const pos=new Map(ids.map((id,i)=>[id,i])),G=ids.map(()=>null),siteId=new Map();
  for(const l of L.slice(1)){const c=splitLine(l,sep),i=pos.get(c[ix[0]]);if(i==null)continue;const a=+c[ix[2]],t=+c[ix[3]];if(!(t>0)||!(a>=0)||a>t)continue;
    const s=c[ix[1]];if(!siteId.has(s))siteId.set(s,siteId.size);(G[i]||(G[i]=new Map())).set(siteId.get(s),[a,t])}
  return G}

// ---------- run ----------
const t0=Date.now();
const table=readTable(pos[0]),meta=readMeta(pos[1]);
for(const c of ['subject','well','plate','household'])if(opt[c]&&!meta.head.includes(opt[c]))die(`column "${opt[c]}" is not in ${path.basename(pos[1])}; columns are: ${meta.head.join(', ')}`);
const depthOf=s=>{let t=0;for(const v of s.val)t+=v;return t};
const minD=opt.type==='counts'?(opt['min-depth']||0):0;
const inMeta=table.samples.filter(s=>meta.rows.has(s.id)&&String(meta.rows.get(s.id)[opt.subject]??'').trim()!=='');
const missing=table.samples.length-inMeta.length;
const shallow=inMeta.filter(s=>depthOf(s)<minD).length;if(shallow){const keep=inMeta.filter(s=>depthOf(s)>=minD);inMeta.length=0;inMeta.push(...keep)}
if(inMeta.length<4)die('fewer than 4 table samples have a subject in the metadata; do the sample IDs match?');
const data={feats:table.feats,samples:inMeta},lab=inMeta.map(s=>String(meta.rows.get(s.id)[opt.subject]).trim());
const col=c=>opt[c]?inMeta.map(s=>meta.rows.get(s.id)[opt[c]]??''):null;
say(`SwapSeq ${VERSION}: ${inMeta.length} samples from ${new Set(lab).size} people${missing?` (${missing} table samples without a subject were left out)`:''}${shallow?` (${shallow} under ${minD} reads left out)`:''}, ${opt.type} mode`);
const M=opt.type==='markers'?markerDistances(data):distances(data);
let res=annotate(data,lab,analyse(data,lab,{prior:opt.prior},M),{house:col('household'),well:col('well'),plate:col('plate')});
const power=opt.noPower?null:measurePower(data,lab,res,M);
let host=null;
if(opt.host){host=hostCheck(readHost(opt.host,inMeta.map(s=>s.id)),lab);res=addHost(res,host);
  say(host.ok?`Host DNA: ${host.used} samples with enough host reads.`:`Host DNA not used: ${host.reason}.`)}
const flags=res.rows.filter(r=>r.level==='flag'),checks=res.rows.filter(r=>r.level==='check');
const pct=v=>v==null||!Number.isFinite(v)?'n/a':Math.round(100*v)+'%';

// ---------- output ----------
const f4=v=>Number.isFinite(v)?(+v).toFixed(4):'';
const cols=['sample','labelled_as','reads','level','kind','probability','looks_like','swap_partner','partner_labelled_as','rare_asvs','rare_shared_with_look_alike','rare_shared_with_own','well','same_household','in_sheet_shift','host_verdict','host_looks_like','host_llr','host_sites','flagged_by_host'];
const lines=[`# SwapSeq ${VERSION} QC audit; table: ${path.basename(pos[0])}; metadata: ${path.basename(pos[1])}; subject column: ${opt.subject}; type: ${opt.type}; date: ${new Date().toISOString().slice(0,10)}`,cols.join('\t')];
for(const r of res.rows){const pair=r.kind==='swap';
  lines.push([r.id,r.subject,opt.type==='counts'?depthOf(data.samples[r.i]):'',r.scorable?r.level:'not_checked',r.level==='ok'?'':r.kind,f4(r.prob),r.looksLike??'',pair?r.partner??'':'',pair?r.partnerSubject??'':'',
    r.rare?r.rare.total:'',r.rare?r.rare.withOther:'',r.rare?r.rare.withOwn:'',r.well??'',r.cohab?'yes':'',r.shift?'yes':'',
    r.host?r.host.verdict:'',r.host?.looksLike??'',r.host&&Number.isFinite(r.host.llr)?r.host.llr.toFixed(1):'',r.host?r.host.sites??'':'',r.byHost?'yes':''].join('\t'))}
fs.writeFileSync(opt.out+'.tsv',lines.join('\n')+'\n');
const methods=`Sample identity was checked with SwapSeq ${VERSION} (${opt.type==='markers'?'strain-marker Jaccard fingerprints':'soft Jaccard and Aitchison fingerprints'}${host&&host.ok?', plus host-DNA genotype matching':''}, calibrated on errors planted in this dataset; expected mislabel rate ${opt.prior*100}%). Of ${res.nScorable} samples from ${res.nSubjects} participants, ${flags.length} ${flags.length===1?'was':'were'} flagged with strong evidence and ${checks.length} listed for checking against laboratory records.${power?` Measured detection in this dataset, from errors hidden in copies of the labels: ${pct(power.swapWithCheck)} of swapped tubes and ${pct(power.singleWithCheck)} of single mislabels (${pct(power.swap)} and ${pct(power.single)} as strong flags).`:''} (Singh 2026, https://doi.org/10.5281/zenodo.23963612)`;
const summary={version:VERSION,type:opt.type,samples:res.rows.length,scorable:res.nScorable,subjects:res.nSubjects,strongFlags:flags.length,checkList:checks.length,
  swapPairs:res.pairs,shifts:res.shifts||[],power,host:host?{ok:host.ok,reason:host.reason||null,samplesWithHostDNA:host.used||0,flaggedByHost:res.rows.filter(r=>r.byHost).length}:null,
  methodsSentence:methods,seconds:(Date.now()-t0)/1000};
fs.writeFileSync(opt.out+'.json',JSON.stringify(summary,null,1));
say(`Strong flags: ${flags.length} · check list: ${checks.length} · swap pairs: ${res.pairs.filter(p=>p.level==='flag').length}`);
if(power)say(`Measured power: ${pct(power.swap)} of hidden swaps caught as strong flags, ${pct(power.swapWithCheck)} counting the check list`);
const low=r=>opt.type==='counts'&&depthOf(data.samples[r.i])<1000;
for(const r of flags.slice(0,20))say(`  FLAG ${r.id} (labelled ${r.subject}) looks like ${r.looksLike??'?'}${r.kind==='swap'&&r.partner?`, swap partner ${r.partner}`:''}${r.byHost?' [host DNA]':''}${low(r)?` [low depth: ${depthOf(data.samples[r.i])} reads]`:''}`);
if(flags.length>20)say(`  ... and ${flags.length-20} more in ${opt.out}.tsv`);
say(`Wrote ${opt.out}.tsv and ${opt.out}.json (${summary.seconds.toFixed(1)} s)`);
process.exit(opt.failOnFlags&&flags.length?2:0);
