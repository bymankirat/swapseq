// Does a plate-aware prior help? Samples are laid out on 96-well plates (random order, or sorted by person), swaps are
// planted between direct neighbours (different people, 3+ samples each) or between random samples, and the pipeline
// is run without and with the plate map.
import * as C from '../../method/core.mjs';import {get} from '../robustness/loadany.mjs';import {load} from '../ye2026_comparison/loadye.mjs';import fs from 'fs';
const key=process.argv[2],REPS=+(process.env.REPS||10),KS=(process.env.KS||'3,5,10').split(',').map(Number);
const ds=key.startsWith('ye:')?load(key.slice(3)):get(key),data=ds.data,lab=ds.lab,n=lab.length,M=C.distances(data);
const cnt=new Map();lab.forEach(l=>cnt.set(l,(cnt.get(l)||0)+1));
function layout(kind,R){let ord=[...Array(n).keys()];if(kind==='random'){for(let i=n-1;i>0;i--){const j=Math.floor(R()*(i+1));[ord[i],ord[j]]=[ord[j],ord[i]]}}else ord.sort((a,b)=>String(lab[a]).localeCompare(String(lab[b]))||a-b);
  const well=new Array(n),plate=new Array(n),pos=new Array(n);ord.forEach((i,k)=>{const p=Math.floor(k/96),w=k%96,c=Math.floor(w/8),r=w%8;well[i]=String.fromCharCode(65+r)+String(c+1).padStart(2,'0');plate[i]='P'+(p+1);pos[i]=[p,r,c]});return {well,plate,pos}}
const NP=n<300?1:Math.max(2,Math.round(n*.003));const out={key,n,res:{}};
const run=(L,lay,k)=>{const r=C.analyse(data,L,{seed:3},M);C.annotate(data,L,r,k?{well:lay.well,plate:lay.plate,plateK:k}:{});return C.promotePairs(r,null)};
for(const kind of ['random','sorted']){for(const mode of ['neighbour','random','none']){
  const acc={};for(const k of [0,...KS])acc[k]={c:0,o:0};let T=0,O=0;
  for(let rep=0;rep<(mode==='none'?3:REPS);rep++){const R=C.rng(7000+rep);for(let q=0;q<10;q++)R();const lay=layout(kind,R),L=lab.slice(),sw=new Set();let g=0;
    const at=new Map();lay.pos.forEach((p,i)=>at.set(p.join(','),i));
    while(mode!=='none'&&sw.size<2*NP&&g++<2e5){const a=Math.floor(R()*n);let b;
      if(mode==='neighbour'){const [p,r,c]=lay.pos[a],d=[[1,0],[-1,0],[0,1],[0,-1]][Math.floor(R()*4)];b=at.get([p,r+d[0],c+d[1]].join(','));if(b==null)continue}else b=Math.floor(R()*n);
      if(a===b||L[a]===L[b]||sw.has(a)||sw.has(b)||cnt.get(lab[a])<3||cnt.get(lab[b])<3)continue;[L[a],L[b]]=[L[b],L[a]];sw.add(a);sw.add(b)}
    T+=sw.size;O+=n-sw.size;
    for(const k of [0,...KS]){const r=run(L,lay,k);for(const x of r.rows)if(x.flag)(sw.has(x.i)?acc[k].c++:acc[k].o++)}}
  out.res[kind+'/'+mode]={T,O,...Object.fromEntries(Object.entries(acc).map(([k,v])=>[k,[v.c,v.o]]))};
  console.log(key.padEnd(18),(kind+'/'+mode).padEnd(17),'hidden',String(T).padStart(3),[0,...KS].map(k=>`k=${k}: ${T?Math.round(100*acc[k].c/T)+'%':'-'} / ${(1000*acc[k].o/O).toFixed(1)}`).join('   '))}}
fs.mkdirSync('out',{recursive:true});fs.writeFileSync(`out/${key.replace(':','_')}.json`,JSON.stringify(out));
