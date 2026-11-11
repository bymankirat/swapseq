// The fair head-to-head again, with the current pipeline (analyse + promotePairs) at several expected error rates.
// Same hidden swaps and seeds as fair.mjs. For each setting: share of hidden swapped samples flagged, and other samples
// flagged per 1,000, for strong flags alone and for strong flags plus the check list.
//   YE_DATA=.../Find-mislabels_Profiles_and_metadata node settings.mjs TABLE [prior,prior,...]
import * as C from '../../method/core.mjs';import {load} from './loadye.mjs';import fs from 'fs';
const name=process.argv[2],priors=(process.argv[3]||'.01,.03,.05,.1').split(',').map(Number),REPS=+(process.env.REPS||10);
const {data,lab}=load(name),n=lab.length,M=C.distances(data);
const cnt=new Map();lab.forEach(l=>cnt.set(l,(cnt.get(l)||0)+1));const pool=[...Array(n).keys()].filter(i=>cnt.get(lab[i])>=3);
const NP=n<300?1:Math.max(3,Math.round(n*.003));const out={name,n,settings:{}};
for(const prior of priors){let c=0,cc=0,o=0,oc=0,T=0,O=0;
  for(let rep=0;rep<REPS;rep++){const R=C.rng(5000+rep),L=lab.slice(),sw=new Set();let g=0;
    while(sw.size<2*NP&&g++<1e5){const a=pool[Math.floor(R()*pool.length)],b=pool[Math.floor(R()*pool.length)];if(a===b||L[a]===L[b]||sw.has(a)||sw.has(b))continue;[L[a],L[b]]=[L[b],L[a]];sw.add(a);sw.add(b)}
    T+=sw.size;O+=n-sw.size;const r=C.promotePairs(C.analyse(data,L,{prior,seed:rep+3},M),null);
    for(const x of r.rows){if(x.flag)(sw.has(x.i)?c++:o++);if(x.level!=='ok')(sw.has(x.i)?cc++:oc++)}}
  out.settings[prior]={hidden:T,others:O,flag:{caught:c,other:o},flagOrCheck:{caught:cc,other:oc}};
  console.log(name.padEnd(16),String(prior).padEnd(5),`flags ${Math.round(100*c/T)}% / ${(1000*o/O).toFixed(1)}   flags + check list ${Math.round(100*cc/T)}% / ${(1000*oc/O).toFixed(1)}`)}
fs.writeFileSync(`settings_${name}.json`,JSON.stringify(out));
