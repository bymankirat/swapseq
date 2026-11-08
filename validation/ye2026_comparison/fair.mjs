// Fair head-to-head: hide swaps among ALL samples of people with 3+ samples (nothing excluded), then count for each
// method the share of hidden swapped samples it flags (recall) and how many OTHER samples it flags (per 1,000).
import fs from 'fs';import {execFileSync} from 'child_process';import * as C from '../../method/core.mjs';import {load} from './loadye.mjs';
const S='.',D=(process.env.YE_DATA||'./Find-mislabels_Profiles_and_metadata')+'/',YE=process.env.YE_TOOL||'./Find-abnormality/find_abnormality.py';
const name=process.argv[2],REPS=+(process.env.REPS||10);
const {data,lab}=load(name),n=lab.length,M=C.distances(data);
const metaRows=fs.readFileSync(D+name+'_meta.txt','utf8').split(/\r?\n/).filter(Boolean);const head=metaRows[0],mrows=metaRows.slice(1).map(l=>l.split('\t'));
const idx=new Map(data.samples.map((s,i)=>[s.id,i]));
const cnt=new Map();lab.forEach(l=>cnt.set(l,(cnt.get(l)||0)+1));const pool=[...Array(n).keys()].filter(i=>cnt.get(lab[i])>=3);
const NP=n<300?1:Math.max(3,Math.round(n*.003));
const wd=`${S}/ye_${name}`;fs.mkdirSync(wd,{recursive:true});
function ye(L,tag){const mp=`${wd}/${tag}_meta.txt`;fs.writeFileSync(mp,head+'\n'+mrows.map(r=>{const i=idx.get(r[0]);return [r[0],i==null?r[1]:L[i],...r.slice(2)].join('\t')}).join('\n')+'\n');
  try{execFileSync('python3',['-I',YE,D+name+'_profiles.txt',mp,'-s',tag],{cwd:wd,stdio:'ignore',timeout:1800000})}catch(e){}
  const s=new Set();const p=`${wd}/${tag}_true_labels.txt`;if(fs.existsSync(p))for(const l of fs.readFileSync(p,'utf8').split('\n').slice(1))if(l.trim())s.add(l.split('\t')[0].split(',')[0]);
  const q=`${wd}/${tag}_cheating_sns.txt`;if(fs.existsSync(q))for(const l of fs.readFileSync(q,'utf8').split('\n'))if(l.includes('=='))l.split('==').forEach(x=>s.add(x.trim()));
  return new Set([...s].map(id=>idx.get(id)).filter(i=>i!=null))}
const M_={ye:{c:0,o:0},s1:{c:0,o:0},s1c:{c:0,o:0},s2:{c:0,o:0},s3:{c:0,o:0}};let T=0,O=0;
for(let rep=0;rep<REPS;rep++){const R=C.rng(5000+rep),L=lab.slice(),sw=new Set();let g=0;
  while(sw.size<2*NP&&g++<1e5){const a=pool[Math.floor(R()*pool.length)],b=pool[Math.floor(R()*pool.length)];if(a===b||L[a]===L[b]||sw.has(a)||sw.has(b))continue;[L[a],L[b]]=[L[b],L[a]];sw.add(a);sw.add(b)}
  T+=sw.size;O+=n-sw.size;
  const add=(k,set)=>{for(const i of set)(sw.has(i)?M_[k].c++:M_[k].o++)};
  add('ye',ye(L,`r${rep}`));
  for(const [k,prior] of [['s1',.01],['s2',.02],['s3',.03]]){const r=C.analyse(data,L,{prior,seed:rep+3},M);add(k,r.rows.filter(x=>x.flag).map(x=>x.i));if(k==='s1')add('s1c',r.rows.filter(x=>x.level!=='ok').map(x=>x.i))}}
const out={name,n,hidden:T,others:O,res:Object.fromEntries(Object.entries(M_).map(([k,v])=>[k,{recall:v.c/T,otherPer1000:1000*v.o/O}]))};
fs.writeFileSync(`${S}/fair_${name}.json`,JSON.stringify(out));
console.log(name.padEnd(16),'n',n,'hidden',T,Object.entries(out.res).map(([k,v])=>`| ${k} ${Math.round(100*v.recall)}% other ${v.otherPer1000.toFixed(1)}`).join(' '));
