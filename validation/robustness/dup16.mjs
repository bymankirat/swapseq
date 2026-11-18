// Duplicate check on 16S count tables: a sample re-sequenced under another person's name (counts jittered, then
// resampled to a new depth; rare taxa drop out by sampling), and planted mislabels that must NOT read as duplicates.
import * as C from '../../method/core.mjs';import {get} from './loadany.mjs';import fs from 'fs';
const key=process.argv[2],REPS=+(process.env.REPS||5),SIG=+(process.env.SIG||0.3);const ds=get(key),data=ds.data,lab=ds.lab,n=lab.length;
function gauss(R){return Math.sqrt(-2*Math.log(R()+1e-12))*Math.cos(2*Math.PI*R())}
function pois(R,l){if(l>30)return Math.max(0,Math.round(l+Math.sqrt(l)*gauss(R)));let L=Math.exp(-l),k=0,p=1;do{k++;p*=R()}while(p>L);return k-1}
function replicate(R,s){const tot=s.val.reduce((a,b)=>a+b,0),f=Math.exp(.5*gauss(R));const I=[],V=[];   // new library: depth x0.4-2.7, each taxon jittered, counts drawn again
  for(let k=0;k<s.idx.length;k++){const v=pois(R,s.val[k]*f*Math.exp(SIG*gauss(R)));if(v>0){I.push(s.idx[k]);V.push(v)}}return {id:s.id+'_dup',idx:Int32Array.from(I),val:Float64Array.from(V)}}
const M0=C.distances(data),pub=C.findDuplicates(M0,lab).pairs,pubS=new Set(pub.map(p=>p.a+'|'+p.b));
const cnt=new Map();lab.forEach(l=>cnt.set(l,(cnt.get(l)||0)+1));const multi=[...Array(n).keys()].filter(i=>cnt.get(lab[i])>=2);
let found=0,tot=0,extra=0,mis=0,misN=0;
for(let rep=0;rep<REPS;rep++){const R=C.rng(900+rep);const samples=data.samples.slice(),L=lab.slice(),pl=[],used=new Set();let g=0;
  while(pl.length<3&&g++<1e4){const a=Math.floor(R()*n),b=Math.floor(R()*n);if(a===b||L[a]===L[b]||used.has(a)||used.has(b))continue;samples[b]=replicate(R,data.samples[a]);used.add(a);used.add(b);pl.push(Math.min(a,b)+'|'+Math.max(a,b))}
  const got=new Set(C.findDuplicates(C.distances({feats:data.feats,samples}),L).pairs.map(p=>p.a+'|'+p.b));
  for(const p of pl){tot++;if(got.has(p))found++}for(const q of got)if(!pl.includes(q)&&!pubS.has(q))extra++;
  const L2=lab.slice(),ms=new Set();g=0;while(multi.length>3&&ms.size<3&&g++<1e4){const x=multi[Math.floor(R()*multi.length)],y=Math.floor(R()*n);if(L2[y]===L2[x]||ms.has(x))continue;L2[x]=L2[y];ms.add(x)}
  if(ms.size){misN+=ms.size;for(const p of C.findDuplicates(M0,L2).pairs)if((ms.has(p.a)||ms.has(p.b))&&!pubS.has(p.a+'|'+p.b))mis++}}
const out={key,n,people:cnt.size,published:pub.length,found,tot,extra,mis,misN};fs.mkdirSync('dup16',{recursive:true});fs.writeFileSync(`dup16/${key.replace(':','_')}.json`,JSON.stringify(out));
console.log(key.padEnd(30),'n',String(n).padStart(5),'| published-label pairs',pub.length,'| planted found',`${found}/${tot}`,'| unplanted new pairs',extra,'| mislabels read as duplicate',`${mis}/${misN}`);
