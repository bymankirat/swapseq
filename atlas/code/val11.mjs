// Validate SwapSeq 1.1 measured power vs the 1.0.1 estimate, on clean labels and on labels with ~1% hidden errors.
import fs from 'fs';import * as L from './load2.mjs';
import {analyse,distances,scorable,rng,measurePower} from '../../swap/core.mjs';
const OUT=(process.env.OUT||'.')+'/val11/';fs.mkdirSync(OUT,{recursive:true});
const all={...L.VALID,...L.NEW,...L.BOREN};for(const n of L.ploegNiches())all['ploeg_'+n.replace(/\W+/g,'_')]=()=>L.ploeg(n);
const key=process.argv[2];let ds=all[key]();
if(!/relative/i.test(ds.type||'')){const keep=[];ds.data.samples.forEach((s,i)=>{if(s.val.reduce((a,b)=>a+b,0)>=1000)keep.push(i)});
  ds={...ds,data:{feats:ds.data.feats,samples:keep.map(i=>ds.data.samples[i])},lab:keep.map(i=>ds.lab[i])}}
const data=ds.data,lab=ds.lab,n=data.samples.length,M=distances(data);
let t=Date.now();const res0=analyse(data,lab,{},M);const tA=(Date.now()-t)/1000;
t=Date.now();const mp0=measurePower(data,lab,res0,M);const tP=(Date.now()-t)/1000;
const truth=JSON.parse(fs.readFileSync('../results/'+key+'.json','utf8')).measured;
// error scenario: K independent evaluation runs, each with its own hidden errors
const K=+(process.env.K||5),ok=scorable(lab),cnt=new Map();lab.forEach(l=>cnt.set(l,(cnt.get(l)||0)+1));
const base=new Set(res0.rows.filter(r=>r.level!=='ok').map(r=>r.i));
const idx=[...Array(n).keys()].filter(i=>ok[i]&&!base.has(i)&&cnt.get(lab[i])>=3);
const small=n<300,NS=small?1:Math.max(3,Math.round(n*.003)),NO=small?1:Math.max(2,Math.round(n*.002));
const E=[];
if(idx.length>=4*(NS+NO))for(let k=0;k<K;k++){const R=rng(50000+k*31),Lb=lab.slice(),used=new Set(),sw=[];
  const pick=()=>{for(;;){const i=idx[Math.floor(R()*idx.length)];if(!used.has(i))return i}};
  while(sw.length<NS){const a=pick(),b=pick();if(a===b||Lb[a]===Lb[b])continue;[Lb[a],Lb[b]]=[Lb[b],Lb[a]];used.add(a);used.add(b);sw.push(a,b)}
  let one=0;while(one<NO){const a=pick(),b=pick();if(a===b||Lb[a]===Lb[b])continue;Lb[a]=Lb[b];used.add(a);used.add(b);one++}
  const r=analyse(data,Lb,{seed:k+11},M);const mp=measurePower(data,Lb,r,M);
  E.push({actualSwap:sw.filter(i=>r.rows[i].flag).length/sw.length,oldSwap:r.power.swap,newSwap:mp?mp.swap:null})}
const mean=a=>a.length?a.reduce((x,y)=>x+y,0)/a.length:null;
const out={key,name:ds.name,n,subjects:new Set(lab).size,analyseSeconds:tA,powerSeconds:tP,reps:mp0&&mp0.reps,
  clean:{truth:truth&&truth.swapCaught,old:res0.power.swap,new:mp0&&mp0.swap,newSingle:mp0&&mp0.single,truthSingle:truth&&truth.relabelCaught},
  withErrors:E.length?{actual:mean(E.map(e=>e.actualSwap)),old:mean(E.map(e=>e.oldSwap)),new:mean(E.filter(e=>e.newSwap!=null).map(e=>e.newSwap)),runs:E.length}:null};
fs.writeFileSync(OUT+key+'.json',JSON.stringify(out,null,1));
const p=v=>v==null||Number.isNaN(v)?' n/a':String(Math.round(v*100)).padStart(4)+'%';
console.log(key.padEnd(30),'n',n,'| clean: truth',p(out.clean.truth),'old',p(out.clean.old),'new',p(out.clean.new),'| 1% errors: actual',p(out.withErrors?.actual),'old',p(out.withErrors?.old),'new',p(out.withErrors?.new),'| +',tP.toFixed(1)+'s');
