// Small studies (validation section 12): subsample K people (up to m samples each) from real public tables, hide ONE swap, compare rules.
import * as C from '../../method/core.mjs';import {get} from './loadall.mjs';import fs from 'fs';
const OLD=await import(process.env.OLD_CORE||'../../method/core.mjs');  // OLD_CORE = path to an earlier core.mjs to compare with
const key=process.argv[2],REPS=+(process.env.REPS||30);const ds=get(key);
const V={old:OLD,new:C};
const by=new Map();ds.lab.forEach((l,i)=>{if(!by.has(l))by.set(l,[]);by.get(l).push(i)});const people=[...by.keys()].filter(p=>by.get(p).length>=3);
const out={key,cfg:{}};
for(const [K,m] of [[6,4],[8,4],[10,4],[15,5],[25,5],[40,6]]){if(people.length<K)continue;const res={};for(const v in V)res[v]={c:0,T:0,o:0,O:0,clean:0,cleanN:0};
  for(let rep=0;rep<REPS;rep++){const R=C.rng(31*K+rep+1);const pp=people.slice().sort(()=>R()-.5).slice(0,K);const idx=[];for(const p of pp){const a=by.get(p).slice().sort(()=>R()-.5).slice(0,m);idx.push(...a)}
    const data={feats:ds.data.feats,samples:idx.map(i=>ds.data.samples[i])},lab=idx.map(i=>ds.lab[i]),n=lab.length,M=C.distances(data);
    let a,b,g=0;do{a=Math.floor(R()*n);b=Math.floor(R()*n)}while((a===b||lab[a]===lab[b])&&g++<1e4);const L=lab.slice();[L[a],L[b]]=[L[b],L[a]];
    for(const v in V){const r0=V[v].promotePairs(V[v].analyse(data,lab,{seed:rep+3},M),null);res[v].clean+=r0.rows.filter(x=>x.flag).length;res[v].cleanN+=n;
      const r=V[v].promotePairs(V[v].analyse(data,L,{seed:rep+3},M),null);res[v].T+=2;res[v].c+=(r.rows[a].flag?1:0)+(r.rows[b].flag?1:0);res[v].O+=n-2;res[v].o+=r.rows.filter((x,i)=>x.flag&&i!==a&&i!==b).length}}
  out.cfg[K+'x'+m]=res;
  console.log(key.padEnd(22),`${K} people x ${m}:`,Object.entries(res).map(([v,x])=>`${v} ${Math.round(100*x.c/x.T)}% (other ${(1000*x.o/x.O).toFixed(1)}, clean ${(1000*x.clean/x.cleanN).toFixed(1)})`).join(' | '))}
fs.mkdirSync('smallfinal',{recursive:true});fs.writeFileSync(`smallfinal/${key}.json`,JSON.stringify(out));
