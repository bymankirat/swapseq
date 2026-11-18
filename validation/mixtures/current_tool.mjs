// What does the tool as it is do with a physical mixture of two people's samples?
import * as C from '../../method/core.mjs';import {get} from '../robustness/loadany.mjs';import fs from 'fs';
const key=process.argv[2],REPS=+(process.env.REPS||8);const ds=get(key);const data=ds.data,lab=ds.lab,n=lab.length;
function gauss(R){return Math.sqrt(-2*Math.log(R()+1e-12))*Math.cos(2*Math.PI*R())}
function pois(R,l){if(l>30)return Math.max(0,Math.round(l+Math.sqrt(l)*gauss(R)));let L=Math.exp(-l),k=0,p=1;do{k++;p*=R()}while(p>L);return k-1}
function mix(R,s,t,a){const Ns=s.val.reduce((x,y)=>x+y,0),Nt=t.val.reduce((x,y)=>x+y,0),m=new Map();
  for(let k=0;k<s.idx.length;k++)m.set(s.idx[k],(1-a)*s.val[k]);for(let k=0;k<t.idx.length;k++)m.set(t.idx[k],(m.get(t.idx[k])||0)+a*t.val[k]*Ns/Nt);
  const I=[],V=[];for(const f of [...m.keys()].sort((x,y)=>x-y)){const v=pois(R,m.get(f));if(v>0){I.push(f);V.push(v)}}return {id:s.id,idx:Int32Array.from(I),val:Float64Array.from(V)}}
const cnt=new Map();lab.forEach(l=>cnt.set(l,(cnt.get(l)||0)+1));const pool=[...Array(n).keys()].filter(i=>cnt.get(lab[i])>=3);
if(pool.length<30){console.log(key,'too few');process.exit(0)}
const NP=n<300?1:Math.max(2,Math.round(n*.003));const out={key,n,res:{}};
for(const a of [0.5,0.3,0.1]){let flag=0,check=0,named=0,T=0,oF=0,O=0;
  for(let rep=0;rep<REPS;rep++){const R=C.rng(6000+rep+Math.round(100*a));for(let k=0;k<20;k++)R();const samples=data.samples.slice(),pl=new Map();let g=0;
    while(pl.size<NP&&g++<1e4){const s=pool[Math.floor(R()*pool.length)],t=pool[Math.floor(R()*pool.length)];if(lab[s]===lab[t]||pl.has(s))continue;samples[s]=mix(R,data.samples[s],data.samples[t],a);pl.set(s,lab[t])}
    const d2={feats:data.feats,samples};const r=C.promotePairs(C.analyse(d2,lab,{seed:rep+3},C.distances(d2)),null);
    for(const x of r.rows){if(pl.has(x.i)){T++;if(x.flag)flag++;if(x.level!=='ok'){check++;if(x.looksLike===pl.get(x.i)||x.partnerSubject===pl.get(x.i))named++}}else{O++;if(x.flag)oF++}}}
  out.res[a]={T,flag,check,named,oF,O}}
fs.mkdirSync('mix3',{recursive:true});fs.writeFileSync(`mix3/${key}.json`,JSON.stringify(out));
console.log(key.padEnd(14),'n',String(n).padStart(5),Object.entries(out.res).map(([a,v])=>`a=${a}: flag ${Math.round(100*v.flag/v.T)}% listed ${Math.round(100*v.check/v.T)}% named ${Math.round(100*v.named/Math.max(1,v.check))}% (${v.T})`).join(' | '));
