#!/usr/bin/env node
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mankirat Singh
// SwapSeq on the command line: the same engine as the web page and the R package, for scripts and pipelines.
// Node.js 18 or newer, no other software. Run `swapseq --help`.
import fs from 'fs';import path from 'path';
import {analyse,annotate,distances,markerDistances,measurePower,promotePairs,findDuplicates,reportHTML,findShifts,hostCheck,addHost} from '../method/core.mjs';

const VERSION='1.1.1';
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
  --out PREFIX        Write PREFIX.tsv (one row per sample), PREFIX.json (summary) and PREFIX.html
                      (a report for your records; open it in a browser and print to PDF). Default: swapseq
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
const read=f=>{if(!fs.existsSync(f))die(`file not found: ${f}`);if(fs.statSync(f).isDirectory())die(`${f} is a folder, not a file`);
  const buf=fs.readFileSync(f);if(buf.length===0)die(`${path.basename(f)} is empty`);
  if(buf[0]===0x89&&buf.slice(1,4).toString('latin1')==='HDF')die(`${path.basename(f)} is a BIOM 2 (HDF5) file. Convert it first: biom convert -i ${path.basename(f)} -o table.tsv --to-tsv`);
  if(buf[0]===0x50&&buf[1]===0x4b)die(`${path.basename(f)} looks like an Excel or zip file. Save the sheet as tab-separated text (.tsv) or .csv first`);
  if(buf.subarray(0,4000).includes(0))die(`${path.basename(f)} is not a text file`);
  return buf.toString('utf8').replace(/^﻿/,'')};
const splitLines=txt=>txt.split(/\r\n|\n|\r/).filter(l=>l.trim()!=='');
// separator: tab, comma or semicolon, whichever splits the header line into most cells
const sepOf=(f,line)=>{const n=c=>splitLine(line,c).length;const t=n('\t'),c=n(','),sc=n(';');return t>=c&&t>=sc?'\t':(sc>c?';':',')};
function splitLine(l,sep){const o=[];let cur='',q=false;for(let i=0;i<l.length;i++){const ch=l[i];if(ch==='"'){if(q&&l[i+1]==='"'){cur+='"';i++}else q=!q;continue}if(ch===sep&&!q){o.push(cur.trim());cur=''}else cur+=ch}o.push(cur.trim());return o}
const uniq=(ids,what)=>{const seen=new Set();for(const id of ids){if(seen.has(id))die(`${what}: "${id}" appears more than once. Every sample needs its own ID`);seen.add(id)}};
const BLANK=new Set(['','na','nan','null','none','n/a','<na>']);
function readTable(f){
  const txt=read(f),name=path.basename(f);
  if(txt.trimStart().startsWith('{')){let j;try{j=JSON.parse(txt)}catch(e){die(`${name} starts like JSON but cannot be read: ${e.message}`)}
    if(!j.rows||!j.columns||!j.data)die('this JSON file is not a BIOM 1.0 table');
    const feats=j.rows.map(r=>r.id),per=j.columns.map(()=>new Map());
    const put=(r,c,v)=>{if(!Number.isFinite(+v)||+v<0)die('the table has negative or non-numeric values');if(+v>0&&per[c])per[c].set(r,(per[c].get(r)||0)+(+v))};
    if(j.matrix_type==='dense')j.data.forEach((row,r)=>row.forEach((v,c)=>put(r,c,v)));else for(const [r,c,v] of j.data)put(r,c,v);
    const ids=j.columns.map(c=>String(c.id).trim());uniq(ids,name);
    return {feats,samples:ids.map((id,c)=>{const a=[...per[c]].sort((x,y)=>x[0]-y[0]);return {id,idx:Int32Array.from(a.map(x=>x[0])),val:Float64Array.from(a.map(x=>x[1]))}})}}
  let L=splitLines(txt);
  L=L.filter(l=>!(l.startsWith('# ')&&!/^#\s*(OTU|ASV|Feature|feature)[ _]?ID/i.test(l)));   // biom "# Constructed from" lines
  if(L.length<2)die(`${name} has a header line but no data rows`);
  const sep=sepOf(f,L[0]);let rows=L.map(l=>splitLine(l,sep));
  // trailing empty cells (a tab at the end of every line) are dropped
  const trail=r=>{let k=r.length;while(k>1&&r[k-1]==='')k--;return k};
  if(rows.every(r=>trail(r)<r.length)){const cut=Math.max(...rows.map(trail));rows=rows.map(r=>r.slice(0,cut))}
  let head=rows[0];const w=rows[1].length;
  // R's write.table() leaves out the corner cell: the header is one cell short of every data row
  if(rows.slice(1).every(r=>r.length===head.length+1)){head=['',...head];rows[0]=head}
  rows.forEach((r,k)=>{if(r.length!==head.length)die(`${name}, line ${k+1} has ${r.length} columns but the header has ${head.length}. Every line needs the same number of columns`)});
  let feats,ids,get;
  if(!opt.rowsAreSamples){feats=rows.slice(1).map(r=>r[0]);ids=head.slice(1);get=(s,fi)=>rows[fi+1][s+1]}
  else{ids=rows.slice(1).map(r=>r[0]);feats=head.slice(1);get=(s,fi)=>rows[s+1][fi+1]}
  // columns (or, with --samples-in-rows, rows) that hold text in every cell, such as a taxonomy column, are not samples
  const isNum=x=>x===''||BLANK.has(x.toLowerCase())||Number.isFinite(+x);
  const textOnly=ids.map((_,s)=>{let num=0,n=0;for(let fi=0;fi<feats.length;fi++){const x=get(s,fi);if(x==='')continue;n++;if(Number.isFinite(+x))num++}return n>0&&num===0});
  const nText=textOnly.filter(Boolean).length;
  if(nText>ids.length/2)die(`${name}: most columns hold text, not numbers. Is it features in rows and samples in columns? (or add --samples-in-rows). A decimal comma (0,5) must be a point (0.5)`);
  const samples=[];
  ids.forEach((id,s)=>{if(textOnly[s])return;const I=[],V=[];for(let fi=0;fi<feats.length;fi++){const raw=get(s,fi);if(BLANK.has(raw.toLowerCase()))continue;const v=+raw;
      if(!Number.isFinite(v))die(`${name}: "${raw}" is not a number (sample ${id}, feature ${feats[fi]})`);if(v<0)die('the table has negative values');if(v>0){I.push(fi);V.push(v)}}
    samples.push({id:String(id),idx:Int32Array.from(I),val:Float64Array.from(V)})});
  uniq(samples.map(x=>x.id),name);
  return {feats,samples,textColumns:nText}}
function readMeta(f){const L=splitLines(read(f)),name=path.basename(f);const sep=sepOf(f,L[0]);
  const head=splitLine(L[0],sep),out=new Map();
  for(const l of L.slice(1)){if(l.startsWith('#'))continue;const c=splitLine(l,sep),o={};head.forEach((h,i)=>o[h]=c[i]??'');
    if(out.has(c[0]))die(`${name}: sample "${c[0]}" appears more than once. Every sample needs one row`);out.set(c[0],o)}
  return {head,rows:out}}
function readHost(f,ids){const L=splitLines(read(f));const sep=sepOf(f,L[0]),h=splitLine(L[0],sep).map(x=>x.toLowerCase());
  const ix=['sample','site','alt','total'].map(k=>h.indexOf(k));if(ix.some(i=>i<0))die('--host file needs the columns sample, site, alt and total');
  const pos=new Map(ids.map((id,i)=>[id,i])),G=ids.map(()=>null),siteId=new Map();
  for(const l of L.slice(1)){const c=splitLine(l,sep),i=pos.get(c[ix[0]]);if(i==null)continue;const a=+c[ix[2]],t=+c[ix[3]];if(!(t>0)||!(a>=0)||a>t)continue;
    const s=c[ix[1]];if(!siteId.has(s))siteId.set(s,siteId.size);(G[i]||(G[i]=new Map())).set(siteId.get(s),[a,t])}
  return G}

// ---------- run ----------
const t0=Date.now();
for(const ext of ['.tsv','.json','.html'])for(const f of [pos[0],pos[1],opt.host].filter(Boolean))if(path.resolve(opt.out+ext)===path.resolve(f))die(`--out ${opt.out} would overwrite the input file ${f}; choose another name`);
const table=readTable(pos[0]),meta=readMeta(pos[1]);
if(table.textColumns)say(`Note: ${table.textColumns} text column${table.textColumns===1?' was':'s were'} not read as samples (taxonomy or notes).`);
for(const c of ['subject','well','plate','household'])if(opt[c]&&!meta.head.includes(opt[c])){const m=meta.head.filter(h=>h.trim().toLowerCase()===String(opt[c]).trim().toLowerCase());if(m.length===1)opt[c]=m[0]}
for(const c of ['subject','well','plate','household'])if(opt[c]&&!meta.head.includes(opt[c]))die(`column "${opt[c]}" is not in ${path.basename(pos[1])}; columns are: ${meta.head.join(', ')}`);
const depthOf=s=>{let t=0;for(const v of s.val)t+=v;return t};
const minD=opt.type==='counts'?(opt['min-depth']||0):0;
const inMeta=table.samples.filter(s=>meta.rows.has(s.id)&&!BLANK.has(String(meta.rows.get(s.id)[opt.subject]??'').trim().toLowerCase()));
const missing=table.samples.length-inMeta.length;
const shallow=inMeta.filter(s=>depthOf(s)<minD).length;if(shallow){const keep=inMeta.filter(s=>depthOf(s)>=minD);inMeta.length=0;inMeta.push(...keep)}
if(inMeta.length<4){const ex=a=>a.slice(0,3).map(x=>`"${x}"`).join(', ');
  die(`fewer than 4 table samples have a subject in the metadata; do the sample IDs match? Table IDs start ${ex(table.samples.map(x=>x.id))}; metadata IDs start ${ex([...meta.rows.keys()])}.${opt.rowsAreSamples?'':' If samples are in rows, add --samples-in-rows.'}`)}
{const c=new Map();for(const s of inMeta){const l=String(meta.rows.get(s.id)[opt.subject]).trim();c.set(l,(c.get(l)||0)+1)}
  if(c.size<2)die(`every sample has the same subject ("${[...c.keys()][0]}"). SwapSeq compares people, so it needs at least 3`);
  if(c.size<3)die('only 2 people in the table. SwapSeq needs at least 3 people to tell a mix-up from a look-alike');
  if(![...c.values()].some(v=>v>=2))die(`no person has two or more samples (column "${opt.subject}" has ${c.size} different values for ${inMeta.length} samples). SwapSeq needs repeat samples from the same people; is this the right column?`)}
const data={feats:table.feats,samples:inMeta},lab=inMeta.map(s=>String(meta.rows.get(s.id)[opt.subject]).trim());
const col=c=>opt[c]?inMeta.map(s=>meta.rows.get(s.id)[opt[c]]??''):null;
say(`SwapSeq ${VERSION}: ${inMeta.length} samples from ${new Set(lab).size} people${missing?` (${missing} table samples without a subject were left out)`:''}${shallow?` (${shallow} under ${minD} reads left out)`:''}, ${opt.type} mode`);
const M=opt.type==='markers'?markerDistances(data):distances(data);
let res=annotate(data,lab,analyse(data,lab,{prior:opt.prior},M),{house:col('household'),well:col('well'),plate:col('plate')});
const power=opt.noPower?null:measurePower(data,lab,res,M);
promotePairs(res,power);
const dups=findDuplicates(M,lab,{},inMeta.map(x=>x.id)).pairs;
let host=null;
if(opt.host){host=hostCheck(readHost(opt.host,inMeta.map(s=>s.id)),lab);res=addHost(res,host);
  say(host.ok?`Host DNA: ${host.used} samples with enough host reads.`:`Host DNA not used: ${host.reason}.`)}
const flags=res.rows.filter(r=>r.level==='flag'),checks=res.rows.filter(r=>r.level==='check');
const pct=v=>v==null||!Number.isFinite(v)?'n/a':Math.round(100*v)+'%';

// ---------- output ----------
const f4=v=>Number.isFinite(v)?(+v).toFixed(4):'';
const cols=['sample','labelled_as','reads','level','kind','probability','looks_like','swap_partner','partner_labelled_as','rare_asvs','rare_shared_with_look_alike','rare_shared_with_own','well','same_household','in_sheet_shift','host_verdict','host_looks_like','host_llr','host_sites','flagged_by_host','flagged_as_pair'];
const lines=[`# SwapSeq ${VERSION} QC audit; table: ${path.basename(pos[0])}; metadata: ${path.basename(pos[1])}; subject column: ${opt.subject}; type: ${opt.type}; date: ${new Date().toISOString().slice(0,10)}`,cols.join('\t')];
for(const r of res.rows){const pair=r.kind==='swap';
  lines.push([r.id,r.subject,opt.type==='counts'?depthOf(data.samples[r.i]):'',(r.scorable||r.byHost)?r.level:'not_checked',r.level==='ok'?'':r.kind,f4(r.prob),r.looksLike??'',pair?r.partner??'':'',pair?r.partnerSubject??'':'',
    r.rare?r.rare.total:'',r.rare?r.rare.withOther:'',r.rare?r.rare.withOwn:'',r.well??'',r.cohab?'yes':'',r.shift?'yes':'',
    r.host?r.host.verdict:'',r.host?.looksLike??'',r.host&&Number.isFinite(r.host.llr)?r.host.llr.toFixed(1):'',r.host?r.host.sites??'':'',r.byHost?'yes':'',r.byPair?'yes':''].join('\t'))}
try{fs.writeFileSync(opt.out+'.tsv',lines.join('\n')+'\n')}catch(e){die(`cannot write ${opt.out}.tsv (${e.code||e.message}). Does the folder exist?`)}
const methods=`Sample identity was checked with SwapSeq ${VERSION} (${opt.type==='markers'?'strain-marker Jaccard fingerprints':'soft Jaccard and Aitchison fingerprints'}${host&&host.ok?', plus host-DNA genotype matching':''}, calibrated on errors planted in this dataset; expected mislabel rate ${opt.prior*100}%). Of ${res.nScorable} samples from ${res.nSubjects} participants, ${flags.length} ${flags.length===1?'was':'were'} flagged with strong evidence and ${checks.length} listed for checking against laboratory records${dups.length?`; ${dups.length} pair${dups.length===1?'':'s'} of samples labelled as different people had near-identical communities (possible duplicates)`:''}.${power?` Measured detection in this dataset, from errors hidden in copies of the labels: ${pct(power.swapWithCheck)} of swapped tubes and ${pct(power.singleWithCheck)} of single mislabels (${pct(power.swap)} and ${pct(power.single)} as strong flags).`:''} (Singh 2026, https://doi.org/10.5281/zenodo.23963612)`;
const manyErrors=res.nScorable&&res.pairs.length>=3&&2*res.pairs.length/res.nScorable>=1.5*opt.prior;
const summary={version:VERSION,type:opt.type,samples:res.rows.length,scorable:res.nScorable,subjects:res.nSubjects,strongFlags:flags.length,checkList:checks.length,
  swapPairs:res.pairs,shifts:res.shifts||[],duplicates:dups.map(p=>({sampleA:p.idA,subjectA:p.subjectA,sampleB:p.idB,subjectB:p.subjectB,distance:+p.distance[0].toFixed(4)})),power,host:host?{ok:host.ok,reason:host.reason||null,samplesWithHostDNA:host.used||0,flaggedByHost:res.rows.filter(r=>r.byHost).length}:null,
  manyErrorsNote:!!manyErrors,methodsSentence:methods,seconds:(Date.now()-t0)/1000};
fs.writeFileSync(opt.out+'.json',JSON.stringify(summary,null,1));
fs.writeFileSync(opt.out+'.html',reportHTML({version:'SwapSeq '+VERSION,file:path.basename(pos[0]),metaFile:path.basename(pos[1]),date:new Date().toISOString().slice(0,10),type:opt.type,prior:opt.prior,
  columns:['subject','well','plate','household'].filter(c=>opt[c]).map(c=>c+'='+opt[c]).join(', '),rows:res.rows,pairs:res.pairs,dups:dups.map(p=>({a:p.idA,sa:p.subjectA,b:p.idB,sb:p.subjectB})),power,nScorable:res.nScorable,nSubjects:res.nSubjects,
  host:host?{ok:host.ok,reason:host.reason,used:host.used}:null,uncalibrated:res.uncalibrated,powerSkipped:!!opt.noPower,methods,notes:[missing?`${missing} table samples without a subject were left out.`:'',shallow?`${shallow} samples under ${minD} reads were left out.`:'',opt.noPower?'Power was not measured (--no-power).':''].filter(Boolean)}));
say(`Strong flags: ${flags.length} · check list: ${checks.length} · swap pairs: ${res.pairs.filter(p=>p.level==='flag').length} · possible duplicates: ${dups.length}`);
if(manyErrors)say(`Note: ${(100*2*res.pairs.length/res.nScorable).toFixed(1)}% of samples sit in swap pairs, more than a ${opt.prior*100}% error rate explains. Try --prior 0.02 or 0.05 (in tests with 3% of samples swapped: ~80% caught instead of 40-55%).`);
if(res.uncalibrated)say('Too few samples with repeats to calibrate: no sample was judged. SwapSeq needs at least 5 samples from people with 2 or more samples.');
else if(!power&&!opt.noPower)say('Power was not measured: it needs people with 3 or more samples. A clean result here does not show the labels are right.');
if(power)say(`Measured power: ${pct(power.swap)} of hidden swaps caught as strong flags, ${pct(power.swapWithCheck)} counting the check list`);
const low=r=>opt.type==='counts'&&depthOf(data.samples[r.i])<1000;
for(const r of flags.slice(0,20))say(`  FLAG ${r.id} (labelled ${r.subject}) looks like ${r.looksLike??'?'}${r.kind==='swap'&&r.partner?`, swap partner ${r.partner}`:''}${r.byHost?' [host DNA]':''}${low(r)?` [low depth: ${depthOf(data.samples[r.i])} reads]`:''}`);
if(flags.length>20)say(`  ... and ${flags.length-20} more in ${opt.out}.tsv`);
for(const d of dups.slice(0,10))say(`  DUPLICATE? ${d.idA} (labelled ${d.subjectA}) and ${d.idB} (labelled ${d.subjectB}) are near-identical: one tube or file under two names, or a sample under the wrong person?`);
say(`Wrote ${opt.out}.tsv, ${opt.out}.json and ${opt.out}.html (${summary.seconds.toFixed(1)} s)`);
process.exit(opt.failOnFlags&&flags.length?2:0);
