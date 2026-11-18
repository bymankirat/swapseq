// old (pushed) core vs new core on every full table: published labels, one swap per ~300 samples, 3% swapped
const OLD=await import(process.env.OLD_CORE||'../../method/core.mjs');import * as NEW from '../../method/core.mjs';import {get} from './loadany.mjs';import fs from 'fs';
const key=process.argv[2],REPS=+(process.env.REPS||6);const d=get(key),data=d.data,lab=d.lab;
const n=lab.length;const cnt=new Map();lab.forEach(l=>cnt.set(l,(cnt.get(l)||0)+1));const pool=[...Array(n).keys()].filter(i=>cnt.get(lab[i])>=3);
const out={key,n};const t0=Date.now();
for(const [nm,C] of [['old',OLD],['new',NEW]]){const M=C.distances(data);const run=(L,seed)=>{const r=C.promotePairs(C.analyse(data,L,{seed},M),null);return r};
  const c=run(lab,7);const o={clean:c.rows.filter(x=>x.flag).length,check:c.rows.filter(x=>x.level!=='ok').length,sep:c.separation,sep2:c.separation2,reid:c.reid,dups:C.findDuplicates(M,lab).pairs.length};
  for(const [s,NS,base] of [['light',n<300?1:Math.max(3,Math.round(n*.003)),5000],['heavy',Math.max(1,Math.round(.03*n/2)),7000]]){let cc=0,T=0,oo=0,O=0,ck=0;
    for(let rep=0;rep<REPS;rep++){const R=C.rng(base+rep),L=lab.slice(),sw=new Set();let g=0;while(sw.size<2*NS&&g++<1e5&&pool.length>4){const a=pool[Math.floor(R()*pool.length)],b=pool[Math.floor(R()*pool.length)];if(a===b||L[a]===L[b]||sw.has(a)||sw.has(b))continue;[L[a],L[b]]=[L[b],L[a]];sw.add(a);sw.add(b)}
      if(!sw.size)break;const r=run(L,rep+3);T+=sw.size;O+=n-sw.size;for(const x of r.rows){if(x.flag)(sw.has(x.i)?cc++:oo++);if(x.level!=='ok'&&!sw.has(x.i))ck++}}
    o[s]={c:cc,T,o:oo,O,ck}}
  out[nm]=o}
fs.mkdirSync('regress',{recursive:true});fs.writeFileSync(`regress/${key.replace(':','_')}.json`,JSON.stringify(out));
const f=(v,s)=>`${Math.round(100*v[s].c/Math.max(1,v[s].T))}%/${(1000*v[s].o/Math.max(1,v[s].O)).toFixed(1)}`;
console.log(key.padEnd(30),'n',String(n).padStart(5),'sep2',out.new.sep2.toFixed(2),'reid',out.new.reid.toFixed(2),'| flags',out.old.clean,'->',out.new.clean,'check',out.old.check,'->',out.new.check,'dups',out.old.dups,'->',out.new.dups,'| light',f(out.old,'light'),'->',f(out.new,'light'),'| heavy',f(out.old,'heavy'),'->',f(out.new,'heavy'),`(${Math.round((Date.now()-t0)/1000)}s)`);
