// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mankirat Singh
// Loads the public validation datasets into {feats, samples} + subject labels.
import fs from 'fs';import zlib from 'zlib';
import * as hdf5 from './lib/hdf5.mjs';
import {fromCSC} from './core.mjs';
const B=(process.env.SWAPSEQ_DATA||'./data')+'/'; // folder holding the cloned source repositories (see README)

export function readTSV(p){const L=fs.readFileSync(p,'utf8').split(/\r?\n/).filter(x=>x&&!x.startsWith('#q2'));const h=L[0].split('\t');return L.slice(1).map(l=>{const c=l.split('\t'),o={};h.forEach((k,i)=>o[k]=c[i]);return o})}
export function biomJSON(txt){const j=JSON.parse(txt);const feats=j.rows.map(r=>r.id),samp=j.columns.map(c=>c.id);
  const per=samp.map(()=>[]);for(const [r,c,v] of j.data)if(v>0)per[c].push([r,v]);
  return {feats,samples:samp.map((id,c)=>{const a=per[c].sort((x,y)=>x[0]-y[0]);return {id,idx:Int32Array.from(a.map(x=>x[0])),val:Float64Array.from(a.map(x=>x[1]))}})}}
export function biomH5(buf){const h=new hdf5.File(buf.buffer.slice(buf.byteOffset,buf.byteOffset+buf.byteLength),'t');
  return fromCSC(h.get('observation/ids').value,h.get('sample/ids').value,h.get('sample/matrix/indptr').value,h.get('sample/matrix/indices').value,h.get('sample/matrix/data').value)}
function subset(d,keep){return {feats:d.feats,samples:d.samples.filter(s=>keep.has(s.id))}}

export function smp(site){
  const d=biomJSON(zlib.gunzipSync(fs.readFileSync(B+`data/smp/otu-tables/otu_table_even10000_${site}.ts_only.biom.gz`)).toString());
  const m=new Map(readTSV(B+'data/smp/analysis-results/issue-13/map.txt').map(r=>[r['#SampleID'],r]));
  const keep=new Set(d.samples.filter(s=>m.has(s.id)&&!/^control/i.test(m.get(s.id).PersonalID)).map(s=>s.id));const x=subset(d,keep);
  return {name:'SMP '+site,data:x,lab:x.samples.map(s=>m.get(s.id).PersonalID),time:x.samples.map(s=>+m.get(s.id).WeeksSinceStart),visit:x.samples.map(s=>m.get(s.id).WeeksSinceStart),meta:m};
}
export function ecam(){
  const f=fs.readdirSync(B+'x/ecam-table-maturity')[0];
  const d=biomH5(fs.readFileSync(B+`x/ecam-table-maturity/${f}/data/feature-table.biom`));
  const m=new Map(readTSV(B+'data/q2-longitudinal/q2_longitudinal/tests/data/ecam_map_maturity.txt').map(r=>[r['#SampleID'],r]));
  const keep=new Set(d.samples.filter(s=>m.has(s.id)).map(s=>s.id));const x=subset(d,keep);
  return {name:'ECAM infants',data:x,lab:x.samples.map(s=>m.get(s.id).studyid),time:x.samples.map(s=>+m.get(s.id).day_of_life),visit:x.samples.map(s=>m.get(s.id).month),meta:m};
}

// generic features-x-samples text table (first column = feature id), optional sample filter
export function tableTSV(p,sep='\t',keepCol=()=>true){
  const L=fs.readFileSync(p,'utf8').replace(/\r/g,'\n').split('\n').filter(x=>x.trim());
  const unq=s=>s.replace(/^"|"$/g,'');const h=L[0].split(sep).map(unq);const cols=[];h.forEach((c,i)=>{if(i>0&&keepCol(c))cols.push(i)});
  const feats=[],per=cols.map(()=>[]);
  L.slice(1).forEach((l,r)=>{const c=l.split(sep);feats.push(unq(c[0]));cols.forEach((ci,k)=>{const v=+c[ci];if(v>0)per[k].push([r,v])})});
  return {feats,samples:cols.map((ci,k)=>({id:h[ci],idx:Int32Array.from(per[k].map(x=>x[0])),val:Float64Array.from(per[k].map(x=>x[1]))}))};
}
export function dietstudy(level='burst'){
  const m=new Map(readTSV(B+'data/dietstudy_analyses/data/maps/SampleID_map.txt').map(r=>[r['#SampleID'],r]));
  const file=level==='burst'?'data/microbiome/BURST.tax.txt':'data/microbiome/processed_sample/taxonomy_counts_s.txt';
  const d=tableTSV(B+'data/dietstudy_analyses/'+file,'\t',c=>m.has(c)&&m.get(c).UserName&&!/^Blank/i.test(m.get(c).UserName)&&m.get(c).SampleType!=='Blank');
  return {name:'Diet study ('+(level==='burst'?'strain-level':'species')+')',data:d,lab:d.samples.map(s=>m.get(s.id).UserName),time:d.samples.map(s=>+m.get(s.id).StudyDayNo),visit:d.samples.map(s=>m.get(s.id).StudyDayNo),meta:m};
}
export function mdsine(cohort){
  const m=new Map(readTSV(B+`data/MDSINE2_Paper/datasets/gibson/${cohort}/raw_tables/metadata.tsv`).map(r=>[r.sampleID,r]));
  const d=tableTSV(B+`data/MDSINE2_Paper/datasets/gibson/${cohort}/raw_tables/counts.tsv`,'\t',c=>m.has(c)&&/^\d+$/.test(m.get(c).subject));
  return {name:`MDSINE2 mice (${cohort} donor)`,data:d,lab:d.samples.map(s=>'mouse'+m.get(s.id).subject),time:d.samples.map(s=>+m.get(s.id).time),visit:d.samples.map(s=>m.get(s.id).time),meta:m};
}
// Immigration Microbiome Project (Vangay et al. 2018): 16S de novo OTUs; every person with 2+ samples
export function imp(){
  const m=new Map(readTSV(B+'data/IMP_analyses/data/map.txt').map(r=>[r['#SampleID'],r]));
  const good=c=>m.has(c)&&!(m.get(c).Exclude||'').trim()&&m.get(c)['Subject.ID'];
  const d0=tableTSV(B+'data/IMP_analyses/data/denovo/final_otu.txt','\t',good);
  const cnt=new Map();d0.samples.forEach(s=>{const k=m.get(s.id)['Subject.ID'];cnt.set(k,(cnt.get(k)||0)+1)});
  const d={feats:d0.feats,samples:d0.samples.filter(s=>cnt.get(m.get(s.id)['Subject.ID'])>=2)};
  return {name:'IMP adults (16S)',data:d,lab:d.samples.map(s=>m.get(s.id)['Subject.ID']),time:d.samples.map(s=>+(m.get(s.id)['Sample.Month']||0)),visit:d.samples.map(s=>m.get(s.id)['Sample.Month']),meta:m};
}
