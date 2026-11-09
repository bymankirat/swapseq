// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mankirat Singh
// SwapSeq Atlas runner: one dataset per call -> results/<key>.json (study-level) + private/<key>.json (sample-level)
import fs from 'fs';import * as L from './load2.mjs';
import {analyse,annotate,distances,scorable,rng,promotePairs,findDuplicates,measurePower} from '../../method/core.mjs';
import {load as loadYe} from '../../validation/ye2026_comparison/loadye.mjs'; // set YE_DATA to the folder of Ye et al. profiles
const REPS_ENV=process.env.REPS?+process.env.REPS:null;
const OUT=process.env.OUT||'./';
const MIN_DEPTH=1000;
const all={...L.VALID,...L.NEW,...L.BOREN};for(const n of L.ploegNiches())all['ploeg_'+n.replace(/\W+/g,'_')]=()=>L.ploeg(n);
all.dietOriginal=()=>{const ds=L.VALID.diet();const lab=ds.lab.slice(),ids=ds.data.samples.map(s=>s.id);const a=ids.indexOf('MCT.f.0086'),b=ids.indexOf('MCT.f.0417');[lab[a],lab[b]]=[lab[b],lab[a]];
  return {...ds,lab,name:'Diet and microbiome, labels before the authors\' correction',control:'known real swap',design:'Same data with one tube swap the authors later corrected'}};
// shotgun tables released by Ye et al. (2026), MetaPhlAn relative abundances (study-level results only)
const ye=(file,name,ref,host,design)=>()=>{const {data,lab}=loadYe(file);return {data,lab,time:lab.map(()=>0),name,ref,type:'Shotgun (MetaPhlAn relative abundance)',host,design,ye:true}};
all.yeBrooks=ye('BrooksB_2017','Preterm infants in NICU (Brooks 2017)','Brooks et al. 2017, Nat Commun; table from Ye et al. 2026','Preterm infant gut','Repeated stool in the first weeks of life');
all.yeT2d=ye('HMP_2019_t2d','Prediabetes, iHMP (Zhou 2019)','Zhou et al. 2019, Nature; table from Ye et al. 2026','Adult gut','Repeated stool over up to four years');
all.yeHall=ye('HallAB_2017','IBD patients and controls (Hall 2017)','Hall et al. 2017, Genome Med; table from Ye et al. 2026','Gut, IBD patients and controls','Repeated stool over about a year');
all.yePRJEB38984=ye('PRJEB38984','Shotgun cohort ENA PRJEB38984','ENA PRJEB38984; table from Ye et al. 2026','Gut (cohort details in Ye et al. 2026)','Repeated stool');
all.yePRJEB70966=ye('PRJEB70966','Shotgun cohort ENA PRJEB70966','ENA PRJEB70966; table from Ye et al. 2026','Gut (cohort details in Ye et al. 2026)','Repeated stool');
all.yeFMT=ye('PRJEB72385_FMT','Faecal transplant cohort ENA PRJEB72385','ENA PRJEB72385; table from Ye et al. 2026','Gut, faecal transplant donors and recipients','Repeated stool around faecal microbiota transplantation');
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
const res=annotate(ds.data,ds.lab,promotePairs(analyse(ds.data,ds.lab,{},M),null));
const dups=findDuplicates(M,ds.lab).pairs.length;
const rp=measurePower(ds.data,ds.lab,res,M);if(rp)promotePairs({rows:[],pairs:[],separation:res.separation},rp);
const flags=res.rows.filter(r=>r.level==='flag'),checks=res.rows.filter(r=>r.level==='check');
// measured detection: hide test errors in the labels (about 1% of samples), rerun, count what is caught (as in the SwapSeq benchmark)
const base0=new Set(res.rows.filter(r=>r.level!=='ok').map(r=>r.i));const cnt=new Map();ds.lab.forEach(l=>cnt.set(l,(cnt.get(l)||0)+1));
// Atlas v0.3 protocol: about 0.5% of samples per repeat (0.3% in swaps, 0.2% one-way relabels). Tables under 500 samples (80 repeats) get a
// single error per repeat (a swap on even repeats, a relabel on odd ones), so no table is tested at more than about 1% errors.
// (v0.2 planted at least 3 swaps and 2 relabels from 300 samples up, up to 2.7% errors in mid-size tables.)
const small=n<500,NSWAP0=small?1:Math.max(1,Math.round(n*0.003)),NONE0=small?1:Math.max(1,Math.round(n*0.002));
const idx=[...Array(n).keys()].filter(i=>ok[i]&&!base0.has(i)&&cnt.get(ds.lab[i])>=3);
let sw=0,sw2=0,swN=0,pairOK=0,ow=0,ow2=0,owN=0,fp=0,fp2=0,fpN=0;
const REPS=REPS_ENV||(n<500?80:n<1500?20:10);
if(idx.length>=4*(NSWAP0+NONE0)&&new Set(idx.map(i=>ds.lab[i])).size>=3)for(let rep=0;rep<REPS;rep++){const R=rng(1000+rep),lab=ds.lab.slice(),used=new Set(),swaps=[],ones=[];
  const NSWAP=small?(rep%2===0?1:0):NSWAP0,NONE=small?(rep%2===1?1:0):NONE0;
  const pick=()=>{for(;;){const i=idx[Math.floor(R()*idx.length)];if(!used.has(i))return i}};
  while(swaps.length<NSWAP){const a=pick(),b=pick();if(a===b||lab[a]===lab[b])continue;[lab[a],lab[b]]=[lab[b],lab[a]];used.add(a);used.add(b);swaps.push([a,b])}
  while(ones.length<NONE){const a=pick(),b=pick();if(a===b||lab[a]===lab[b])continue;lab[a]=lab[b];used.add(a);used.add(b);ones.push(a)}
  const r=promotePairs(analyse(ds.data,lab,{seed:rep+3},M),null);
  for(const [a,b] of swaps){for(const x of [a,b]){swN++;if(r.rows[x].flag)sw++;if(r.rows[x].level!=='ok')sw2++}
    if(r.pairs.some(p=>(p.a===r.rows[a].id&&p.b===r.rows[b].id)||(p.a===r.rows[b].id&&p.b===r.rows[a].id)))pairOK++}
  for(const a of ones){owN++;if(r.rows[a].flag)ow++;if(r.rows[a].level!=='ok')ow2++}
  const planted=new Set(swaps.flat().concat(ones));
  for(let i=0;i<n;i++){if(planted.has(i)||base0.has(i)||!r.rows[i].scorable)continue;fpN++;if(r.rows[i].flag)fp++;if(r.rows[i].level!=='ok')fp2++}}
const measured=swN?{reps:REPS,protocol:'v0.3',swapsPerRep:small?0.5:NSWAP0,relabelsPerRep:small?0.5:NONE0,swapCaught:sw/swN,swapCaughtWithReview:sw2/swN,swapPairNamed:pairOK/(swN/2),relabelCaught:owN?ow/owN:null,relabelCaughtWithReview:owN?ow2/owN:null,falseFlagsPer1000:1000*fp/fpN,falseChecksPer1000:1000*fp2/fpN}:null;
const per=new Map();ds.lab.forEach(l=>per.set(l,(per.get(l)||0)+1));const spp=[...per.values()].filter(v=>v>=2).sort((a,b)=>a-b);
const out={key,name:ds.name,ref:ds.ref,type:ds.type,host:ds.host,design:ds.design,validation:!!ds.val,control:ds.control||null,
  samples:n,droppedLowDepth:nLow,subjects:per.size,scorable:res.nScorable,subjectsWith2:spp.length,medianSamplesPerSubject:spp[spp.length>>1]||0,
  reidentified:reN?reHit/reN:null,selfEstimatedPowerOnOriginal:res.power,measured,
  separation:res.separation,promotedOnOriginal:res.promoted||0,possibleDuplicates:dups,engine:'1.1.1 November update',reportedPower:rp?{swap:rp.swap,swapBeforePromotion:rp.swapBeforePromotion??rp.swap,swapWithCheck:rp.swapWithCheck}:null,strongFlags:flags.length,checkList:checks.length,strongSwapPairs:res.pairs.filter(p=>p.level==='flag').length,allSwapPairs:res.pairs.length,
  flagsPer1000:res.nScorable?1000*flags.length/res.nScorable:null,checksPer1000:res.nScorable?1000*checks.length/res.nScorable:null,
  householdOrRareNote:flags.filter(r=>r.rare&&r.rare.withOther>r.rare.withOwn).length,
  seconds:(Date.now()-t0)/1000};
fs.mkdirSync(OUT+'results',{recursive:true});fs.mkdirSync(OUT+'private',{recursive:true});
fs.writeFileSync(OUT+'results/'+key+'.json',JSON.stringify(out,null,1));
fs.writeFileSync(OUT+'private/'+key+'.json',JSON.stringify(res.rows.filter(r=>r.level!=='ok').map(r=>({id:r.id,subject:r.subject,level:r.level,kind:r.kind,prob:r.prob,looksLike:r.looksLike,partner:r.partner,partnerSubject:r.partnerSubject,time:ds.time[r.i],rare:r.rare||null})),null,1));
console.log(key,n,'samples',out.subjects,'subj','reID',(out.reidentified*100).toFixed(1)+'%','caught',measured?(measured.swapCaught*100).toFixed(0)+'%':'n/a','review',measured?(measured.swapCaughtWithReview*100).toFixed(0)+'%':'','relabel',measured&&measured.relabelCaught!=null?(measured.relabelCaught*100).toFixed(0)+'%':'','fp/1000',measured?measured.falseFlagsPer1000.toFixed(2):'','flags',flags.length,'checks',checks.length,'pairs',res.pairs.length,out.seconds.toFixed(1)+'s');
