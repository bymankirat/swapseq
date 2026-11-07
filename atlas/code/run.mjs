// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mankirat Singh
// SwapSeq Atlas runner: one dataset per call -> results/<key>.json (study-level) + private/<key>.json (sample-level)
import fs from 'fs';import * as L from './load2.mjs';
import {analyse,annotate,distances,scorable,rng} from '../../swap/core.mjs';
const REPS=+(process.env.REPS||10);
const OUT='/tmp/claude-0/-home-claude/d674b933-5c36-5fd7-9c63-4fe2a591bf95/scratchpad/atlas/';
const MIN_DEPTH=1000;
const all={...L.VALID,...L.NEW,...L.BOREN};for(const n of L.ploegNiches())all['ploeg_'+n.replace(/\W+/g,'_')]=()=>L.ploeg(n);
all.dietOriginal=()=>{const ds=L.VALID.diet();const lab=ds.lab.slice(),ids=ds.data.samples.map(s=>s.id);const a=ids.indexOf('MCT.f.0086'),b=ids.indexOf('MCT.f.0417');[lab[a],lab[b]]=[lab[b],lab[a]];
  return {...ds,lab,name:'Diet and microbiome, labels before the authors\' correction',control:'known real swap',design:'Same data with one tube swap the authors later corrected'}};
const key=process.argv[2];if(key==='list'){console.log(Object.keys(all).join(' '));process.exit(0)}
let ds=all[key]();const t0=Date.now();
// uniform QC: drop samples below MIN_DEPTH reads (count data only), then drop subjects left with one sample? (kept: they are unscorable but still act as look-alike candidates)
const relative=/relative/i.test(ds.type||'');let nLow=0;
if(!relative){const keep=[];ds.data.samples.forEach((s,i)=>{const d=s.val.reduce((a,b)=>a+b,0);if(d>=MIN_DEPTH)keep.push(i);else nLow++});
  ds={...ds,data:{feats:ds.data.feats,samples:keep.map(i=>ds.data.samples[i])},lab:keep.map(i=>ds.lab[i]),time:keep.map(i=>ds.time[i])}}
const n=ds.data.samples.length,M=distances(ds.data);
// re-identification: is the most similar other sample from the same subject? (fingerprint distance, scorable samples)
const ok=scorable(ds.lab);let reHit=0,reN=0;
for(let i=0;i<n;i++){if(!ok[i])continue;let best=Infinity,bj=-1;for(let j=0;j<n;j++){if(j===i)continue;const d=M.J[i*n+j];if(d<best){best=d;bj=j}}reN++;if(ds.lab[bj]===ds.lab[i])reHit++}
const res=annotate(ds.data,ds.lab,analyse(ds.data,ds.lab,{},M));
const flags=res.rows.filter(r=>r.level==='flag'),checks=res.rows.filter(r=>r.level==='check');
// measured detection: hide test errors in the labels (about 1% of samples), rerun, count what is caught (as in the SwapSeq benchmark)
const base0=new Set(res.rows.filter(r=>r.level!=='ok').map(r=>r.i));const cnt=new Map();ds.lab.forEach(l=>cnt.set(l,(cnt.get(l)||0)+1));
const small=n<300,NSWAP=small?1:Math.max(3,Math.round(n*0.003)),NONE=small?1:Math.max(2,Math.round(n*0.002));
const idx=[...Array(n).keys()].filter(i=>ok[i]&&!base0.has(i)&&cnt.get(ds.lab[i])>=3);
let sw=0,sw2=0,swN=0,pairOK=0,ow=0,ow2=0,owN=0,fp=0,fp2=0,fpN=0;
if(idx.length>=4*(NSWAP+NONE)&&new Set(idx.map(i=>ds.lab[i])).size>=3)for(let rep=0;rep<REPS;rep++){const R=rng(1000+rep),lab=ds.lab.slice(),used=new Set(),swaps=[],ones=[];
  const pick=()=>{for(;;){const i=idx[Math.floor(R()*idx.length)];if(!used.has(i))return i}};
  while(swaps.length<NSWAP){const a=pick(),b=pick();if(a===b||lab[a]===lab[b])continue;[lab[a],lab[b]]=[lab[b],lab[a]];used.add(a);used.add(b);swaps.push([a,b])}
  while(ones.length<NONE){const a=pick(),b=pick();if(a===b||lab[a]===lab[b])continue;lab[a]=lab[b];used.add(a);used.add(b);ones.push(a)}
  const r=analyse(ds.data,lab,{seed:rep+3},M);
  for(const [a,b] of swaps){for(const x of [a,b]){swN++;if(r.rows[x].flag)sw++;if(r.rows[x].level!=='ok')sw2++}
    if(r.pairs.some(p=>(p.a===r.rows[a].id&&p.b===r.rows[b].id)||(p.a===r.rows[b].id&&p.b===r.rows[a].id)))pairOK++}
  for(const a of ones){owN++;if(r.rows[a].flag)ow++;if(r.rows[a].level!=='ok')ow2++}
  const planted=new Set(swaps.flat().concat(ones));
  for(let i=0;i<n;i++){if(planted.has(i)||base0.has(i)||!r.rows[i].scorable)continue;fpN++;if(r.rows[i].flag)fp++;if(r.rows[i].level!=='ok')fp2++}}
const measured=swN?{reps:REPS,swapsPerRep:NSWAP,relabelsPerRep:NONE,swapCaught:sw/swN,swapCaughtWithReview:sw2/swN,swapPairNamed:pairOK/(swN/2),relabelCaught:owN?ow/owN:null,relabelCaughtWithReview:owN?ow2/owN:null,falseFlagsPer1000:1000*fp/fpN,falseChecksPer1000:1000*fp2/fpN}:null;
const per=new Map();ds.lab.forEach(l=>per.set(l,(per.get(l)||0)+1));const spp=[...per.values()].filter(v=>v>=2).sort((a,b)=>a-b);
const out={key,name:ds.name,ref:ds.ref,type:ds.type,host:ds.host,design:ds.design,validation:!!ds.val,control:ds.control||null,
  samples:n,droppedLowDepth:nLow,subjects:per.size,scorable:res.nScorable,subjectsWith2:spp.length,medianSamplesPerSubject:spp[spp.length>>1]||0,
  reidentified:reN?reHit/reN:null,selfEstimatedPowerOnOriginal:res.power,measured,
  strongFlags:flags.length,checkList:checks.length,strongSwapPairs:res.pairs.filter(p=>p.level==='flag').length,allSwapPairs:res.pairs.length,
  flagsPer1000:res.nScorable?1000*flags.length/res.nScorable:null,checksPer1000:res.nScorable?1000*checks.length/res.nScorable:null,
  householdOrRareNote:flags.filter(r=>r.rare&&r.rare.withOther>r.rare.withOwn).length,
  seconds:(Date.now()-t0)/1000};
fs.mkdirSync(OUT+'results',{recursive:true});fs.mkdirSync(OUT+'private',{recursive:true});
fs.writeFileSync(OUT+'results/'+key+'.json',JSON.stringify(out,null,1));
fs.writeFileSync(OUT+'private/'+key+'.json',JSON.stringify(res.rows.filter(r=>r.level!=='ok').map(r=>({id:r.id,subject:r.subject,level:r.level,kind:r.kind,prob:r.prob,looksLike:r.looksLike,partner:r.partner,partnerSubject:r.partnerSubject,time:ds.time[r.i],rare:r.rare||null})),null,1));
console.log(key,n,'samples',out.subjects,'subj','reID',(out.reidentified*100).toFixed(1)+'%','caught',measured?(measured.swapCaught*100).toFixed(0)+'%':'n/a','review',measured?(measured.swapCaughtWithReview*100).toFixed(0)+'%':'','relabel',measured&&measured.relabelCaught!=null?(measured.relabelCaught*100).toFixed(0)+'%':'','fp/1000',measured?measured.falseFlagsPer1000.toFixed(2):'','flags',flags.length,'checks',checks.length,'pairs',res.pairs.length,out.seconds.toFixed(1)+'s');
