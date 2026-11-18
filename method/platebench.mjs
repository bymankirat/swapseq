// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mankirat Singh
// Plate-pattern errors on the diet study's real plating order (Run1.Plate + SampleOrder):
//  neighbour swap = the labels of two tubes next to each other on a plate are exchanged
//  shift run      = a sample sheet slips by one row: L consecutive tubes each get the next tube's label
import {analyse,annotate,distances,rng} from './core.mjs';import {dietstudy} from './load.mjs';import fs from 'fs';
const ds=dietstudy(process.argv[2]||'burst'),S=ds.data.samples,n=S.length,REPS=+(process.argv[3]||30),L=+(process.argv[4]||5);
const pos=S.map(s=>{const r=ds.meta.get(s.id);return r['Run1.Plate']!=='NA'&&r.SampleOrder!=='NA'?{plate:r['Run1.Plate'],order:+r.SampleOrder}:null});
const byPlate=new Map();pos.forEach((p,i)=>{if(!p)return;if(!byPlate.has(p.plate))byPlate.set(p.plate,[]);byPlate.get(p.plate).push(i)});
for(const a of byPlate.values())a.sort((x,y)=>pos[x].order-pos[y].order);
const runs=[...byPlate.values()];const M=distances(ds.data);
const base0=new Set(analyse(ds.data,ds.lab,{},M).rows.filter(r=>r.level!=='ok').map(r=>r.i));
const res={};
for(const mode of ['neighbour swap','shift run']){let shiftHit=0,falseShift=0,hit=0,hitC=0,N=0,runHit=0,runN=0,fp=0,fpc=0,fpN=0;
  for(let rep=0;rep<REPS;rep++){const R=rng(500+rep),lab=ds.lab.slice(),touched=new Set();
    if(mode==='neighbour swap'){let k=0;while(k<3){const pl=runs[Math.floor(R()*runs.length)],j=Math.floor(R()*(pl.length-1)),a=pl.at(j),b=pl.at(j+1);
        if(lab[a]===lab[b]||touched.has(a)||touched.has(b)||base0.has(a)||base0.has(b))continue;[lab[a],lab[b]]=[lab[b],lab[a]];touched.add(a);touched.add(b);k++}}
    else{for(;;){const pl=runs[Math.floor(R()*runs.length)],j=Math.floor(R()*(pl.length-L));if(j<0)continue;const seg=pl.slice(j,j+L+1);if(seg.some(i=>base0.has(i)))continue;
        for(let q=0;q<L;q++){const a=seg.at(q),b=seg.at(q+1);lab[a]=ds.lab[b];if(lab[a]!==ds.lab[a])touched.add(a)}break}}
    const r=annotate(ds.data,lab,analyse(ds.data,lab,{seed:rep+11},M),{well:S.map((s,i)=>pos[i]?String(pos[i].order):''),plate:S.map((s,i)=>pos[i]?pos[i].plate:'')});let rh=0;
    const tset=[...touched].map(i=>S[i].id);const detected=r.shifts.filter(x=>x.ids.filter(id=>tset.includes(id)).length>=Math.min(3,x.ids.length)).length>0;shiftHit+=(mode==='shift run'&&detected)?1:0;falseShift+=r.shifts.filter(x=>!x.ids.some(id=>tset.includes(id))).length;
    for(const i of touched){N++;if(r.rows[i].level==='flag')hit++;if(r.rows[i].level!=='ok'){hitC++;rh++}}
    if(mode==='shift run'){runN++;if(rh>=Math.ceil(touched.size/2))runHit++}
    for(let i=0;i<n;i++){if(touched.has(i)||base0.has(i)||!r.rows[i].scorable)continue;fpN++;if(r.rows[i].level==='flag')fp++;if(r.rows[i].level!=='ok')fpc++}}
  res[mode]={samplesFlagged:+(hit/N).toFixed(3),withCheck:+(hitC/N).toFixed(3),...(mode==='shift run'?{runsWithHalfCaught:+(runHit/runN).toFixed(3),shiftNamed:+(shiftHit/runN).toFixed(3),runLength:L}:{}),falseShiftReportsPerRepeat:+(falseShift/REPS).toFixed(3),falseFlagsPer1000:+(1000*fp/fpN).toFixed(2),falseChecksPer1000:+(1000*fpc/fpN).toFixed(2),errorsPlanted:N}}
console.log(ds.name,'plates',runs.length,'tubes with position',runs.flat().length,JSON.stringify(res,null,1));
fs.writeFileSync(`plate_results_${process.argv[2]||'burst'}.json`,JSON.stringify({dataset:ds.name,res},null,1));
// where was the real swap on the plate?
const ia=S.findIndex(s=>s.id==='MCT.f.0086'),ib=S.findIndex(s=>s.id==='MCT.f.0417');console.log('real swap positions',JSON.stringify(pos[ia]),JSON.stringify(pos[ib]));
