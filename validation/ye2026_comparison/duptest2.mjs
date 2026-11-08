// Duplicate check variants: planted technical duplicates found, and planted MISLABELS wrongly called duplicates.
import fs from 'fs';import * as C from '../../method/core.mjs';import {load} from './loadye.mjs';
const names=process.argv.slice(2),REPS=+(process.env.REPS||5),SIG=+(process.env.SIG||0.3);
const V=[['global',{}],['pair .5',{same:'pair',f:.5}],['pair .7',{same:'pair',f:.7}],['pair .35',{same:'pair',f:.35}]];
function gauss(R){return Math.sqrt(-2*Math.log(R()+1e-12))*Math.cos(2*Math.PI*R())}
function replicate(R,s){const I=[],Vv=[];const tot=s.val.reduce((a,b)=>a+b,0);
  for(let k=0;k<s.idx.length;k++){const v=s.val[k];if(v/tot<1e-3&&R()<.3)continue;I.push(s.idx[k]);Vv.push(Math.max(1,Math.round(v*Math.exp(SIG*gauss(R)))))}
  return {id:s.id+'_dup',idx:Int32Array.from(I),val:Float64Array.from(Vv)}}
const tot={};
for(const name of names){const {data,lab}=load(name),n=lab.length;if(n<20)continue;const M0=C.distances(data);
  const out={};for(const [k,o] of V){out[k]={pub:C.findDuplicates(M0,lab,o).pairs.length,found:0,tot:0,extra:0,mis:0,misN:0}}
  const cnt=new Map();lab.forEach(l=>cnt.set(l,(cnt.get(l)||0)+1));const multi=[...Array(n).keys()].filter(i=>cnt.get(lab[i])>=2);
  for(let rep=0;rep<REPS;rep++){const R=C.rng(900+rep);const samples=data.samples.slice(),L=lab.slice(),pl=[],used=new Set();
    while(pl.length<3){const a=Math.floor(R()*n),b=Math.floor(R()*n);if(a===b||L[a]===L[b]||used.has(a)||used.has(b))continue;samples[b]=replicate(R,data.samples[a]);used.add(a);used.add(b);pl.push(Math.min(a,b)+'|'+Math.max(a,b))}
    // mislabels: 3 samples of people with repeats relabelled to another person (no duplicate exists)
    const L2=lab.slice(),mis=new Set();let g=0;while(multi.length>3&&mis.size<3&&g++<1e4){const x=multi[Math.floor(R()*multi.length)],y=Math.floor(R()*n);if(L2[y]===L2[x]||mis.has(x))continue;L2[x]=L2[y];mis.add(x)}
    const M=C.distances({feats:data.feats,samples}),pubS=new Map(V.map(([k,o])=>[k,new Set(C.findDuplicates(M0,lab,o).pairs.map(p=>p.a+'|'+p.b))]));
    for(const [k,o] of V){const got=new Set(C.findDuplicates(M,L,o).pairs.map(p=>p.a+'|'+p.b));
      for(const p of pl){out[k].tot++;if(got.has(p))out[k].found++}
      for(const q of got)if(!pl.includes(q)&&!pubS.get(k).has(q))out[k].extra++;
      if(mis.size){const r2=C.findDuplicates(M0,L2,o).pairs;out[k].misN+=mis.size;for(const p of r2)if((mis.has(p.a)||mis.has(p.b))&&!pubS.get(k).has(p.a+'|'+p.b))out[k].mis++}}}
  console.log(name.padEnd(16),'n',n,'|',V.map(([k])=>`${k}: pub ${out[k].pub} found ${out[k].found}/${out[k].tot} extra ${out[k].extra} mislabel-as-dup ${out[k].mis}/${out[k].misN}`).join(' || '));
  for(const [k] of V){tot[k]=tot[k]||{};for(const f in out[k])tot[k][f]=(tot[k][f]||0)+out[k][f]}}
console.log('TOTAL',JSON.stringify(tot));
