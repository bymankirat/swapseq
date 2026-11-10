// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mankirat Singh
// More public tables for the SwapSeq Atlas (version 0.4): Knights lab MLRepo (github.com/knights-lab/MLRepo),
// closed-reference Greengenes OTU tables with QIIME mapping files. Set MLREPO to the cloned repository.
import fs from 'fs';
const ML=(process.env.MLREPO||'./MLRepo')+'/datasets/';
function readMap(p){const L=fs.readFileSync(p,'utf8').split(/\r?\n/).filter(x=>x.trim()),h=L[0].split('\t');return L.slice(1).map(l=>{const c=l.split('\t'),o={};h.forEach((k,i)=>o[k]=c[i]);return o})}
// taxa in rows, samples in columns; only the wanted sample columns are kept (the HMP table has 6,407 columns)
function otuTable(p,want){const txt=fs.readFileSync(p,'utf8');let pos=txt.indexOf('\n');const head=txt.slice(0,pos).replace(/\r$/,'').split('\t');
  const cols=[];head.forEach((id,j)=>{if(j>0&&want.has(id))cols.push(j)});const per=cols.map(()=>({I:[],V:[]})),feats=[];let fi=0;
  while(pos<txt.length){const nx=txt.indexOf('\n',pos+1),line=txt.slice(pos+1,nx<0?txt.length:nx).replace(/\r$/,'');pos=nx<0?txt.length:nx;if(!line)continue;
    const c=line.split('\t');let any=false;for(let k=0;k<cols.length;k++){const v=+c[cols[k]];if(v>0){per[k].I.push(fi);per[k].V.push(v);any=true}}
    if(any){feats.push(c[0]);fi++}else{/* feature absent from the kept samples: index not used */}}
  return {feats,samples:cols.map((j,k)=>({id:head[j],idx:Int32Array.from(per[k].I),val:Float64Array.from(per[k].V)}))}}
function ml(dir,keep,subj,info){const map=readMap(ML+dir+'/mapping-orig.txt').filter(keep),by=new Map(map.map(r=>[r['#SampleID'],r]));
  const d=otuTable(ML+dir+'/gg/otutable.txt',new Set(by.keys()));
  return {data:d,lab:d.samples.map(s=>by.get(s.id)[subj]),time:d.samples.map(()=>0),house:info.house?d.samples.map(s=>by.get(s.id)[info.house]):null,...info}}
const hmp=(site,label,host)=>()=>ml('hmp',r=>r.SPECIFIC_BODY_SITE===site&&r.MISLABELED!=='TRUE','HOST_SUBJECT_ID',
  {name:`HMP1 ${label} (HMP Consortium 2012)`,ref:'Human Microbiome Project Consortium 2012, Nature; table from MLRepo',type:'16S (V3-V5, closed-reference OTUs)',host,design:'Healthy adults, one to three visits over up to two years'});
export const MORE={
  hmpStool:hmp('UBERON:feces','stool','Adult gut'),
  hmpTongue:hmp('UBERON:mucosa of tongue','tongue','Oral: tongue'),
  hmpPlaque:hmp('UBERON:gingiva','dental plaque','Oral: dental plaque'),
  hmpNostril:hmp('UBERON:nostril','nostril','Nose: anterior nares'),
  hmpSkin:hmp('UBERON:skin','skin behind the ear','Skin: behind the ear'),
  turnbaugh:()=>ml('turnbaugh',r=>true,'HOST_SUBJECT_ID',{name:'Obese and lean twins (Turnbaugh 2009)',ref:'Turnbaugh et al. 2009, Nature; table from MLRepo',type:'16S (V2, closed-reference OTUs)',host:'Adult gut, twins and mothers',design:'Twin pairs and their mothers, two stool samples about two months apart',house:'FAMILY'}),
  cafe:()=>ml('bushman_cafe',r=>true,'HOST_SUBJECT_ID',{name:'Controlled feeding, COMBO-CAFE (Wu 2011)',ref:'Wu et al. 2011, Science; table from MLRepo',type:'16S (V1-V2, closed-reference OTUs)',host:'Adult gut',design:'Ten inpatients on a fixed diet, stool daily for ten days'}),
};
