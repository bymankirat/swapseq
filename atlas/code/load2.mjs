// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mankirat Singh
// Loaders for the SwapSeq Atlas: public longitudinal microbiome datasets -> {name, data:{feats,samples}, lab, time, meta}
import fs from 'fs';
import {readTSV,biomJSON,tableTSV,smp,ecam,dietstudy,mdsine,imp} from '../../method/load.mjs';
const A=(process.env.ATLAS_SRC||'./src')+'/'; // folder holding the cloned source repositories
const csvSplit=l=>{const o=[];let cur='',q=false;for(let i=0;i<l.length;i++){const ch=l[i];if(ch==='"'){q=!q;continue}if(ch===','&&!q){o.push(cur);cur=''}else cur+=ch}o.push(cur);return o};
function readCSV(p,skip=0){const L=fs.readFileSync(p,'utf8').replace(/^﻿/,'').split(/\r?\n/).slice(skip).filter(x=>x.trim());const h=csvSplit(L[0]);return L.slice(1).map(l=>{const c=csvSplit(l),o={};h.forEach((k,i)=>o[k]=c[i]);return o})}
// samples-in-rows table (first column sample id, then features) -> {feats,samples}
function rowsTable(p,sep,keep=()=>true,scale=1,skipCols=1){
  const L=fs.readFileSync(p,'utf8').split(/\r?\n/).filter(x=>x.trim());const unq=s=>s.replace(/^"|"$/g,'');
  const h=L[0].split(sep).map(unq);const feats=h.slice(skipCols);const samples=[];
  for(const l of L.slice(1)){const c=l.split(sep);const id=unq(c[0]);if(!keep(id,c))continue;const I=[],V=[];
    for(let k=skipCols;k<c.length;k++){let v=+c[k];if(scale!==1)v=Math.round(v*scale);if(v>0){I.push(k-skipCols);V.push(v)}}
    samples.push({id,idx:Int32Array.from(I),val:Float64Array.from(V)})}
  return {feats,samples}}
function scaleCols(d,scale){for(const s of d.samples){const I=[],V=[];for(let k=0;k<s.idx.length;k++){const v=Math.round(s.val[k]*scale);if(v>0){I.push(s.idx[k]);V.push(v)}}s.idx=Int32Array.from(I);s.val=Float64Array.from(V)}return d}
function dropEmpty(ds){const keep=ds.data.samples.map((s,i)=>s.idx.length>0?i:-1).filter(i=>i>=0);
  return {...ds,data:{feats:ds.data.feats,samples:keep.map(i=>ds.data.samples[i])},lab:keep.map(i=>ds.lab[i]),time:keep.map(i=>ds.time[i])}}
const mk=(name,d,lab,time,extra={})=>dropEmpty({name,data:d,lab,time,...extra});

// ---- MITRE example data (Bogart et al. 2019): samples in rows, counts
function mitre(dir){const m=new Map(fs.readFileSync(A+`gerberlab_mitre/mitre/example_data/${dir}/sample_metadata.csv`,'utf8').split(/\r?\n/).filter(Boolean).map(l=>{const c=l.split(',');return [c[0],{subject:c[1],time:+c[2]}]}));
  const d=rowsTable(A+`gerberlab_mitre/mitre/example_data/${dir}/abundance.csv`,',',id=>m.has(id));return {m,d}}
export function karelia(){const {m,d}=mitre('karelia');return mk('DIABIMMUNE infants (Vatanen 2016)',d,d.samples.map(s=>m.get(s.id).subject),d.samples.map(s=>m.get(s.id).time),{ref:'Vatanen et al. 2016, Cell',type:'16S',host:'Infant gut',design:'3 countries, monthly stool, birth to 3 years'})}
export function digiulio(){const {m,d}=mitre('digiulio');return mk('Pregnancy vaginal (DiGiulio 2015)',d,d.samples.map(s=>m.get(s.id).subject),d.samples.map(s=>m.get(s.id).time),{ref:'DiGiulio et al. 2015, PNAS',type:'16S',host:'Vaginal',design:'Weekly swabs through pregnancy'})}
export function david(){const {m,d}=mitre('david');return mk('Diet switch (David 2014)',d,d.samples.map(s=>m.get(s.id).subject),d.samples.map(s=>m.get(s.id).time),{ref:'David et al. 2014, Nature',type:'16S',host:'Adult gut',design:'Daily stool before, during and after 5-day plant or animal diet'})}

// ---- parafac4microbiome (van der Ploeg): Shao 2019 Baby Biome, van der Ploeg 2024 gingivitis
export function shao(){const P=A+'GRvanderPloeg_parafac4microbiome/data-raw/Shao2019/';
  const meta=new Map(readCSV(P+'Shao2019_sampleMetadata.csv',2).map(r=>[r.Accession,r]));
  const run2s=new Map(readTSV(P+'Shao2019_mapping.tsv').map(r=>[r.run_accession,r.secondary_sample_accession]));
  const L=fs.readFileSync(P+'Shao2019_counts.csv','utf8').split(/\r?\n/).filter(Boolean);const h=L[0].split(',').slice(1).map(x=>x.replace('_taxonomy_profile.txt',''));
  const cols=[];h.forEach((run,k)=>{const ers=run2s.get(run);if(ers&&meta.has(ers))cols.push([k,ers])});
  const per=cols.map(()=>[]);const feats=[];
  L.slice(1).forEach((l,r)=>{const c=l.split(',');feats.push(c[0]||('f'+r));cols.forEach(([k],j)=>{const v=+c[k+1];if(v>0)per[j].push([r,v])})});
  const d={feats,samples:cols.map(([k,ers],j)=>({id:ers,idx:Int32Array.from(per[j].map(x=>x[0])),val:Float64Array.from(per[j].map(x=>x[1]))}))};
  return mk('Baby Biome infants (Shao 2019)',d,d.samples.map(s=>meta.get(s.id).Individual),d.samples.map(s=>+meta.get(s.id).Time_point),{ref:'Shao et al. 2019, Nature',type:'Shotgun (mOTU counts)',host:'Infant gut',design:'Neonatal stool days 4-21 and infancy'})}
export function ploeg(niche){const P=A+'GRvanderPloeg_parafac4microbiome/data-raw/vanderPloeg2024/20221005_wp2/count-table.tsv';
  const L=fs.readFileSync(P,'utf8').split(/\r?\n/).filter(Boolean);const h=L[0].split('\t');const feats=h.slice(5);const samples=[],lab=[],time=[];
  for(const l of L.slice(1)){const c=l.split('\t');if(c[4]!==niche)continue;const I=[],V=[];for(let k=5;k<c.length;k++){const v=+c[k];if(v>0){I.push(k-5);V.push(v)}}
    samples.push({id:c[0],idx:Int32Array.from(I),val:Float64Array.from(V)});lab.push(c[1]);time.push(+c[2])}
  return mk(`Gingivitis challenge, ${niche} (van der Ploeg 2024)`,{feats,samples},lab,time,{ref:'van der Ploeg et al. 2024',type:'16S',host:'Oral: '+niche,design:'Experimental gingivitis, 7 visits over 6 weeks'})}
export function ploegNiches(){const P=A+'GRvanderPloeg_parafac4microbiome/data-raw/vanderPloeg2024/20221005_wp2/count-table.tsv';
  const s=new Set();for(const l of fs.readFileSync(P,'utf8').split(/\r?\n/).slice(1))if(l)s.add(l.split('\t',5)[4]);return [...s]}

// ---- GPZ longitudinal collection
export function taur(){const P=A+'gpz/Taur_2018/';const m=new Map();
  for(const r of readCSV(P+'PRJNA491657_metadata.csv'))if(r.library_strategy==='AMPLICON'&&r.patientID)m.set(r.run_alias.replace('.fastq','')+'..'+r.pool,r);
  const d=tableTSV(P+'otu_counts.tsv','\t',c=>m.has(c));
  return mk('Stem-cell transplant patients (Taur 2018)',d,d.samples.map(s=>'p'+m.get(s.id).patientID),d.samples.map(s=>+m.get(s.id).day),{ref:'Taur et al. 2018, Sci Transl Med',type:'16S',host:'Adult gut (hospital, antibiotics)',design:'Stool from 25 days before to 100 days after allo-HSCT'})}
export function hagan(){const P=A+'gpz/Hagan_2019/';const m=new Map(readCSV(P+'PRJNA505336_metadata.csv').map(r=>[r.run_accession,r]));
  const L=fs.readFileSync(P+'rel-table-ASV_with-DADA2-tax.tsv','utf8').split(/\r?\n/).filter(Boolean);const unq=s=>s.replace(/^"|"$/g,'');const h=L[0].split('\t').map(unq);
  const cols=[];h.forEach((c,i)=>{if(m.has(c))cols.push(i)});const feats=[],per=cols.map(()=>[]);
  L.slice(1).forEach((l,r)=>{const c=l.split('\t');feats.push(unq(c[0]));cols.forEach((ci,k)=>{const v=Math.round(+c[ci]*1e5);if(v>0)per[k].push([r,v])})});
  const d={feats,samples:cols.map((ci,k)=>({id:h[ci],idx:Int32Array.from(per[k].map(x=>x[0])),val:Float64Array.from(per[k].map(x=>x[1]))}))};
  return mk('Antibiotics and vaccines (Hagan 2019)',d,d.samples.map(s=>m.get(s.id).isolate),d.samples.map(s=>+m.get(s.id)['Day of study']),{ref:'Hagan et al. 2019, Cell',type:'16S (relative abundance)',host:'Adult gut',design:'Antibiotics before flu vaccination, stool to day 180'})}
export function palleja(){const P=A+'gpz/Palleja_2018/';const d0=tableTSV(P+'annotated.mOTU.rel_abund.tsv','\t',c=>/^ERAS\d+_Dag/.test(c));const d=scaleCols(d0,1e5);
  return mk('Antibiotic recovery (Palleja 2018)',d,d.samples.map(s=>s.id.split('_')[0]),d.samples.map(s=>+(s.id.match(/Dag(\d+)/)||[0,0])[1]),{ref:'Palleja et al. 2018, Nat Microbiol',type:'Shotgun (mOTU relative)',host:'Adult gut',design:'4-day antibiotic cocktail, stool to day 180'})}

// ---- MicrobeDS: Song et al. 2016 storage study (Qiita 10394)
export function storage(){const P=A+'twbattaglia_MicrobeDS/data-raw/qa10394/';
  const m=new Map(readTSV(P+'10394_prep_1360_qiime_20161216-085956.txt').filter(r=>r.sample_type==='stool').map(r=>[r['#SampleID'],r]));
  const d0=biomJSON(fs.readFileSync(P+'otu_table_json.biom','utf8'));const d={feats:d0.feats,samples:d0.samples.filter(s=>m.has(s.id))};
  return mk('Stool storage conditions (Song 2016)',d,d.samples.map(s=>m.get(s.id).host_subject_id),d.samples.map(s=>0),{ref:'Song et al. 2016, mSystems',type:'16S',host:'Human and dog stool',design:'Aliquots of one stool per donor, many preservation methods, up to 8 weeks'})}
// ---- Mu et al. mouse VRE (lecaolab)
export function muVRE(){const P=A+'lecaolab_Longitudinal_microbiome_data_analysis/Clustering/Case Study/Data/';
  const L=fs.readFileSync(P+'Mu_etal_mouseVRE_FEATURE_TABLE_case_corrected.tsv','utf8').split(/\r?\n/).filter(l=>l&&!l.startsWith('# Constructed'));
  fs.writeFileSync((process.env.TMPDIR||'/tmp')+'/swapseq_mu.tsv',L.join('\n'));
  const d=tableTSV((process.env.TMPDIR||'/tmp')+'/swapseq_mu.tsv','\t',c=>/^11737\.D\d\dM[A-Z]\d$/.test(c));
  return mk('Mouse VRE colonisation (Mu 2021)',d,d.samples.map(s=>s.id.slice(9)),d.samples.map(s=>+s.id.slice(7,9)),{ref:'lecaolab case study (mouse VRE)',type:'16S (batch-corrected)',host:'Mouse gut',design:'Mice sampled daily, days 0-5'})}

// ---- the five validation studies (SwapSeq paper)
const tag=(ds,o)=>({...ds,...o});
export const VALID={
  smpGut:()=>tag(smp('gut'),{name:'Student Microbiome gut',ref:'Flores et al. 2014, Genome Biol',type:'16S',host:'Adult gut',design:'Weekly self-sampling for 10 weeks',val:true}),
  smpTongue:()=>tag(smp('tongue'),{name:'Student Microbiome tongue',ref:'Flores et al. 2014, Genome Biol',type:'16S',host:'Oral: tongue',design:'Weekly self-sampling for 10 weeks',val:true}),
  smpPalm:()=>tag(smp('palm'),{name:'Student Microbiome palm',ref:'Flores et al. 2014, Genome Biol',type:'16S',host:'Skin: palm',design:'Weekly self-sampling for 10 weeks',val:true}),
  smpForehead:()=>tag(smp('forehead'),{name:'Student Microbiome forehead',ref:'Flores et al. 2014, Genome Biol',type:'16S',host:'Skin: forehead',design:'Weekly self-sampling for 10 weeks',val:true}),
  ecam:()=>tag(ecam(),{name:'ECAM infants (Bokulich 2016)',ref:'Bokulich et al. 2016, Sci Transl Med',type:'16S',host:'Infant gut',design:'Monthly stool, birth to 2 years',val:true}),
  diet:()=>tag(dietstudy('burst'),{name:'Diet and microbiome (Johnson 2019)',ref:'Johnson et al. 2019, Cell Host Microbe',type:'Shotgun (strain-level)',host:'Adult gut',design:'Daily stool for 17 days',val:true}),
  mdsH:()=>tag(mdsine('healthy'),{name:'Gnotobiotic mice, healthy donor (Gibson 2025)',ref:'Gibson et al. 2025 (MDSINE2)',type:'16S',host:'Mouse gut',design:'Mice colonised with one donor community',val:true}),
  mdsU:()=>tag(mdsine('uc'),{name:'Gnotobiotic mice, UC donor (Gibson 2025)',ref:'Gibson et al. 2025 (MDSINE2)',type:'16S',host:'Mouse gut',design:'Mice colonised with one donor community',val:true}),
  imp:()=>tag(imp(),{name:'Immigration Microbiome Project (Vangay 2018)',ref:'Vangay et al. 2018, Cell',type:'16S',host:'Adult gut',design:'Immigrants sampled before and after moving to the US',val:true}),
};
export const NEW={karelia,digiulio,david,shao,taur,hagan,palleja,storage,muVRE};

// ---- Borenstein lab curated gut microbiome-metabolome collection (Muller et al. 2022): uniform reprocessing
function boren(dir,name,info,timeCol){const P=A+`borenstein/data/processed_data/${dir}/`;
  const meta=new Map(readTSV(P+'metadata.tsv').map(r=>[r.Sample,r]));
  const f=fs.existsSync(P+'species.counts.tsv')?'species.counts.tsv':'genera.counts.tsv';
  const d=rowsTable(P+f,'\t',id=>meta.has(id)&&meta.get(id).Subject);
  return mk(name,d,d.samples.map(s=>meta.get(s.id).Subject),d.samples.map(s=>{const v=meta.get(s.id)[timeCol];return isNaN(+v)?v:+v}),{...info,level:f.startsWith('species')?'species':'genus'})}
export const BOREN={
  heInfants:()=>boren('HE_INFANTS_MFGM_2019','Infant formula trial (He 2019)',{ref:'He et al. 2019, Sci Rep',type:'16S (genus counts)',host:'Infant gut',design:'Stool at 2, 4, 6 and 12 months'},'Age'),
  kostic:()=>boren('KOSTIC_INFANTS_DIABETES_2015','Infants at risk of type 1 diabetes (Kostic 2015)',{ref:'Kostic et al. 2015, Cell Host Microbe',type:'16S (genus counts)',host:'Infant gut',design:'Monthly stool, birth to 3 years (subset with metabolomics)'},'Age'),
  mars:()=>boren('MARS_IBS_2020','Irritable bowel syndrome (Mars 2020)',{ref:'Mars et al. 2020, Cell',type:'Shotgun (species counts)',host:'Adult gut',design:'Stool over 6 months, including flares'},'Timepoint'),
  poyet:()=>boren('POYET_BIO_ML_2019','BIO-ML stool donors (Poyet 2019)',{ref:'Poyet et al. 2019, Nat Med',type:'Shotgun (genus counts)',host:'Adult gut',design:'Healthy donors sampled repeatedly over up to 2 years'},'Collection_Date'),
  wandro:()=>boren('WANDRO_PRETERMS_2018','Preterm infants in NICU (Wandro 2018)',{ref:'Wandro et al. 2018, mSphere',type:'16S (genus counts)',host:'Preterm infant gut',design:'Stool over the first weeks in intensive care'},'day'),
  ibdmdb:()=>boren('iHMP_IBDMDB_2019','IBD multi-omics, HMP2 (Lloyd-Price 2019)',{ref:'Lloyd-Price et al. 2019, Nature',type:'Shotgun (species counts)',host:'Adult and child gut',design:'Stool every 2 weeks for 1 year (subset with metabolomics)'},'week_num'),
};
