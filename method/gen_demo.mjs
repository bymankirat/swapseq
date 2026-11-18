// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mankirat Singh
// Demo: 16 students from the Student Microbiome Project (gut, weekly), with 2 planted errors.
import {smp} from './load.mjs';import {analyse,rng} from './core.mjs';import fs from 'fs';
const ds=smp('gut');const by=new Map();ds.lab.forEach((l,i)=>{if(!by.has(l))by.set(l,[]);by.get(l).push(i)});
const people=[...by.entries()].filter(([l,a])=>a.length>=9).map(([l])=>l).sort().slice(0,16);
const keepI=[];for(const p of people)keepI.push(...by.get(p).sort((a,b)=>ds.time[a]-ds.time[b]));
const S=keepI.map(i=>ds.data.samples[i]);
const prev=new Map();for(const s of S)for(const f of s.idx)prev.set(f,(prev.get(f)||0)+1);
const feats=[...prev.entries()].filter(([f,c])=>c>=2).map(([f])=>f).sort((a,b)=>a-b),pos=new Map(feats.map((f,k)=>[f,k]));
const samples=S.map(s=>{const I=[],V=[];s.idx.forEach((f,k)=>{if(pos.has(f)){I.push(pos.get(f));V.push(s.val[k])}});return {id:s.id,idx:I,val:V}});
let lab=keepI.map(i=>ds.lab[i]);const time=keepI.map(i=>ds.time[i]);
// rename people so the demo carries no study IDs
const alias=new Map(people.map((p,k)=>[p,'P'+String(k+1).padStart(2,'0')]));lab=lab.map(l=>alias.get(l));
// plant: one tube swap between two people at the same week, and one sample labelled as the wrong person
const R=rng(42);const L=lab.slice();let planted=null;
for(let t=0;t<500&&!planted;t++){const a=Math.floor(R()*L.length),b=Math.floor(R()*L.length);if(L[a]===L[b]||time[a]!==time[b])continue;
  const c=Math.floor(R()*L.length),d=Math.floor(R()*L.length);if(new Set([a,b,c]).size<3||L[c]===L[d]||[a,b].includes(d))continue;
  const L2=L.slice();[L2[a],L2[b]]=[L2[b],L2[a]];L2[c]=L[d];
  const r=analyse({feats,samples:samples.map(s=>({id:s.id,idx:Int32Array.from(s.idx),val:Float64Array.from(s.val)}))},L2);
  const fl=r.rows.filter(x=>x.flag).map(x=>x.i).sort();if(fl.length===3&&fl.includes(a)&&fl.includes(b)&&fl.includes(c)&&r.pairs.length===1){planted={a,b,c,L2,r}}}
const {a,b,c,L2,r}=planted;
console.log('swap',samples[a].id,lab[a],'<->',samples[b].id,lab[b],'week',time[a],' relabel',samples[c].id,lab[c],'->',L2[c],' power',r.power);
const meta='sample-id\tparticipant\tweek\tbody_site\n'+samples.map((s,i)=>`${s.id}\t${L2[i]}\t${time[i]}\tgut`).join('\n')+'\n';
// compact sparse table: one line per sample "id<TAB>delta-encoded idx:count ..."
const tab=samples.map(s=>{let p=0;return s.id+'\t'+s.idx.map((x,k)=>{const d=x-p;p=x;return d.toString(36)+':'+s.val[k].toString(36)}).join(' ')}).join('\n');
fs.writeFileSync('demo.json',JSON.stringify({nFeat:feats.length,table:tab,meta,answer:{swap:[samples[a].id,samples[b].id],relabel:samples[c].id}}));
console.log('features',feats.length,'samples',samples.length,'bytes',fs.statSync('demo.json').size);
