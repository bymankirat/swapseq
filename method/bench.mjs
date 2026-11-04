// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mankirat Singh
// End-to-end benchmark of the shipped pipeline: plant real-looking errors into public longitudinal 16S studies
// and count what SwapSeq reports at its default settings (prior 1%, flag at probability >= 0.5).
import {analyse,distances,rng,scorable} from './core.mjs';
import {smp,ecam,dietstudy,mdsine,imp} from './load.mjs';
import fs from 'fs';
const REPS=+(process.argv[2]||20),PRIOR=+(process.argv[3]||0.01),BASE=process.argv[4]==='bc',CAL=process.env.CAL||'logit',VIS=process.env.VISIT==='1';
function twoPerSubject(ds){ // paired design: keep the first and last sample of each person (like before/after studies)
  const by=new Map();ds.data.samples.forEach((s,i)=>{const l=ds.lab[i];if(!by.has(l))by.set(l,[]);by.get(l).push(i)});
  const keep=[];for(const [l,a] of by){if(a.length<2)continue;a.sort((x,y)=>ds.time[x]-ds.time[y]);keep.push(a[0],a[a.length-1])}
  keep.sort((a,b)=>a-b);return {name:ds.name+' (2 per person)',data:{feats:ds.data.feats,samples:keep.map(i=>ds.data.samples[i])},lab:keep.map(i=>ds.lab[i]),time:keep.map(i=>ds.time[i]),visit:keep.map(i=>ds.visit[i])};
}
const MK={gut:()=>smp('gut'),tongue:()=>smp('tongue'),palm:()=>smp('palm'),forehead:()=>smp('forehead'),ecam,diet:()=>dietstudy('burst'),dietsp:()=>dietstudy('species'),mdh:()=>mdsine('healthy'),mdu:()=>mdsine('uc'),imp};
const want=(process.env.SETS||'gut,tongue,palm,forehead,ecam,gut2,forehead2').split(',');
const sets=want.map(k=>k.endsWith('2')&&MK[k.slice(0,-1)]?twoPerSubject(MK[k.slice(0,-1)]()):MK[k]());
const out=[];
for(const ds of sets){
  const n=ds.data.samples.length,M=distances(ds.data,()=>{},BASE),OPT={calib:CAL,...(VIS?{visit:ds.visit}:{})},ok=scorable(ds.lab);
  // baseline flags on the untouched data: these are excluded from the false-flag count
  const r00=analyse(ds.data,ds.lab,{prior:PRIOR,...OPT},M),base0=new Set(r00.rows.filter(r=>r.flag).map(r=>r.i)),checks0=r00.rows.filter(r=>r.check).length;
  const cnt=new Map();ds.lab.forEach(l=>cnt.set(l,(cnt.get(l)||0)+1));
  const two=ds.name.includes('2 per');const small=!two&&n<300,NSWAP=two||small?1:3,NONE=two?0:small?1:2; // small tables get 1 swap + 1 relabel // about 1% of samples carry an error
  const idx=[...Array(n).keys()].filter(i=>ok[i]&&!base0.has(i)&&(two||cnt.get(ds.lab[i])>=3));
  let sw2=0,ow2=0,fp2=0,sw=0,swN=0,pairOK=0,ow=0,owN=0,fp=0,fpN=0,pw={single:0,swap:0},subjHit=0,subjN=0,subjFP=0,subjFPN=0;
  for(let rep=0;rep<REPS;rep++){const R=rng(1000+rep),lab=ds.lab.slice(),used=new Set(),swaps=[],ones=[];
    const pick=()=>{for(;;){const i=idx[Math.floor(R()*idx.length)];if(!used.has(i))return i}};
    while(swaps.length<NSWAP){const a=pick(),b=pick();if(a===b||lab[a]===lab[b])continue;[lab[a],lab[b]]=[lab[b],lab[a]];used.add(a);used.add(b);swaps.push([a,b])}
    while(ones.length<NONE){const a=pick(),b=pick();if(a===b||lab[a]===lab[b])continue;lab[a]=lab[b];used.add(a);used.add(b);ones.push(a)} // b is just a donor label; mark used so it stays clean
    const r=analyse(ds.data,lab,{seed:rep+3,prior:PRIOR,...OPT},M);pw.single+=r.power.single/REPS;pw.swap+=r.power.swap/REPS;
    for(const [a,b] of swaps){for(const x of [a,b]){swN++;if(r.rows[x].flag)sw++;if(r.rows[x].level!=='ok')sw2++}
      if(r.pairs.some(p=>(p.a===r.rows[a].id&&p.b===r.rows[b].id)||(p.a===r.rows[b].id&&p.b===r.rows[a].id)))pairOK++}
    for(const a of ones){owN++;if(r.rows[a].flag)ow++;if(r.rows[a].level!=='ok')ow2++}
    const planted=new Set(swaps.flat().concat(ones));
    const invP=new Set(two?swaps.flat().map(i=>ds.lab[i]):[]); // with 2 per person, the swapped people's other sample is affected too
    for(let i=0;i<n;i++){if(planted.has(i)||base0.has(i)||!r.rows[i].scorable||invP.has(ds.lab[i]))continue;fpN++;if(r.rows[i].flag)fp++;if(r.rows[i].level!=='ok')fp2++}
    if(two){ // with 2 samples per person a swap makes both people inconsistent; score at the level of people
      const involved=new Set(swaps.flat().map(i=>ds.lab[i]));const flaggedSubj=new Set();
      r.rows.forEach((x,i)=>{if(x.flag){flaggedSubj.add(ds.lab[i])}});
      for(const s of involved){subjN++;if(flaggedSubj.has(s))subjHit++}
      for(const s of new Set(ds.lab)){if(involved.has(s))continue;subjFPN++;if(flaggedSubj.has(s))subjFP++}}
  }
  const res={dataset:ds.name,samples:n,people:new Set(ds.lab).size,swappedSamplesCaught:+(sw/swN).toFixed(3),swapPairsNamedExactly:+(pairOK/(swN/2)).toFixed(3),
    oneWayCaught:owN?+(ow/owN).toFixed(3):null,falseFlagsPer1000:+(1000*fp/fpN).toFixed(2),withReview:{swaps:+(sw2/swN).toFixed(3),oneWay:owN?+(ow2/owN).toFixed(3):null,falsePer1000:+(1000*fp2/fpN).toFixed(2)},selfEstimatedPower:{single:+pw.single.toFixed(3),swap:+pw.swap.toFixed(3)},
    ...(two?{peopleInSwapFlagged:+(subjHit/subjN).toFixed(3),uninvolvedPeopleFlagged:+(subjFP/subjFPN).toFixed(3)}:{}),flagsOnOriginal:base0.size,checksOnOriginal:checks0,reps:REPS,prior:PRIOR};
  out.push(res);console.log(JSON.stringify(res));
}
fs.writeFileSync(`bench_results_prior${PRIOR}${BASE?'_braycurtis':''}${process.env.TAG||''}.json`,JSON.stringify(out,null,1));
