// Experiment: can a sample that is a physical mixture of two people be told from the table alone?
// Statistic per sample s (label A): for every other person B, the share of B's personal taxa (core in B, core in few
// people, never seen in A's other samples) that turn up in s. A mixture shows one person B far above the rest.
import * as C from '../../method/core.mjs';import {get} from '../robustness/loadany.mjs';
const key=process.argv[2],REPS=+(process.env.REPS||3),MINC=+(process.env.MINC||2),NPL=+(process.env.NPL||40);
const ds=get(key);let data=ds.data;const lab=ds.lab,n=lab.length;
if(/relative/i.test(ds.type||''))data={feats:data.feats,samples:C.asCounts(data.samples).S};
function gauss(R){return Math.sqrt(-2*Math.log(R()+1e-12))*Math.cos(2*Math.PI*R())}
function pois(R,l){if(l>30)return Math.max(0,Math.round(l+Math.sqrt(l)*gauss(R)));let L=Math.exp(-l),k=0,p=1;do{k++;p*=R()}while(p>L);return k-1}
function mix(R,s,t,a){const Ns=s.val.reduce((x,y)=>x+y,0),Nt=t.val.reduce((x,y)=>x+y,0),m=new Map();
  for(let k=0;k<s.idx.length;k++)m.set(s.idx[k],(1-a)*s.val[k]);for(let k=0;k<t.idx.length;k++)m.set(t.idx[k],(m.get(t.idx[k])||0)+a*t.val[k]*Ns/Nt);
  const I=[],V=[];for(const f of [...m.keys()].sort((x,y)=>x-y)){const v=pois(R,m.get(f));if(v>0){I.push(f);V.push(v)}}return {id:s.id,idx:Int32Array.from(I),val:Float64Array.from(V)}}
const pres=s=>{const o=[];for(let k=0;k<s.idx.length;k++)if(s.val[k]>=MINC)o.push(s.idx[k]);return o};
// people
const people=[...new Set(lab)],pi=new Map(people.map((p,i)=>[p,i])),P=people.length,L=lab.map(l=>pi.get(l)),np=new Int32Array(P);L.forEach(p=>np[p]++);
const cnt=Array.from({length:P},()=>new Map());const PR=data.samples.map(pres);
PR.forEach((fs,i)=>{const m=cnt[L[i]];for(const f of fs)m.set(f,(m.get(f)||0)+1)});
const donors=[...Array(P).keys()].filter(p=>np[p]>=2);
const isCore=(p,f,minus=0)=>{const c=(cnt[p].get(f)||0)-minus,m=np[p]-(minus?1:0);return c>=2&&c>=0.5*m};
const coreOf=new Map();for(const p of donors)for(const [f,c] of cnt[p])if(c>=2&&c>=0.5*np[p]){if(!coreOf.has(f))coreOf.set(f,[]);coreOf.get(f).push(p)}
const maxR=Math.max(2,Math.floor(0.1*donors.length));
const rareCore=Array.from({length:P},()=>[]);for(const [f,ps] of coreOf)if(ps.length<=maxR)for(const p of ps)rareCore[p].push(f);
// score one sample given its presence list, its label person A, whether it is one of A's counted samples, and optionally a donor sample to hold out
function score(fs,A,minus,holdP=-1,holdSet=null){
  const own=cnt[A],fset=new Set(fs);const inOwn=f=>{const c=(own.get(f)||0)-(minus.has(f)?1:0);return c>=1};
  const Cs=[],who=[];
  for(const B of donors){if(B===A)continue;let U=0,h=0;
    for(const f of rareCore[B]){if(inOwn(f))continue;if(B===holdP&&holdSet.has(f)&&!isCore(B,f,1))continue;U++;if(fset.has(f))h++}
    if(U>=8){Cs.push(h/U);who.push(B)}}
  if(Cs.length<5)return null;const srt=Cs.slice().sort((a,b)=>a-b),med=srt[srt.length>>1],mad=srt.map(x=>Math.abs(x-med)).sort((a,b)=>a-b)[srt.length>>1];
  let bi=0;for(let i=1;i<Cs.length;i++)if(Cs[i]>Cs[bi])bi=i;const second=srt[srt.length-2];
  return {best:Cs[bi],who:who[bi],diff:Cs[bi]-med,z:(Cs[bi]-med)/(1.4826*mad+0.02),gap:Cs[bi]-second,med}}
const obs=[];for(let i=0;i<n;i++){if(np[L[i]]<3)continue;const r=score(PR[i],L[i],new Set(PR[i]));if(r)obs.push({i,...r})}
const stats=['best','diff','z','gap'];const out={key,n,P,scored:obs.length,maxR,rare:rareCore.reduce((a,b)=>a+b.length,0)/Math.max(1,donors.length)};
if(obs.length<30){console.log(key,'too few scorable',obs.length);process.exit(0)}
const q=(arr,p)=>{const s=arr.slice().sort((a,b)=>a-b);return s[Math.min(s.length-1,Math.floor(p*s.length))]};
const thr={};for(const st of stats){const v=obs.map(o=>o[st]);thr[st]={max:Math.max(...v),q99:q(v,.99),q3:v.slice().sort((a,b)=>b-a)[Math.min(v.length-1,3)]}}
const elig=obs.map(o=>o.i),res={};
for(const a of [0.5,0.3,0.2,0.1]){const hit={};for(const st of stats)hit[st]={max:0,q3:0,q99:0,named:0};let tot=0;
  for(let rep=0;rep<REPS;rep++){const R=C.rng(4000+rep*17+Math.round(a*100));for(let k=0;k<20;k++)R();
    for(let k=0;k<NPL;k++){const s=elig[Math.floor(R()*elig.length)];let t,g=0;do{t=Math.floor(R()*n)}while((L[t]===L[s]||np[L[t]]<3)&&g++<1e4);if(L[t]===L[s]||np[L[t]]<3)continue;
      const m=mix(R,data.samples[s],data.samples[t],a),r=score(pres(m),L[s],new Set(PR[s]),L[t],new Set(PR[t]));if(!r)continue;tot++;
      for(const st of stats){if(r[st]>thr[st].max){hit[st].max++;if(r.who===L[t])hit[st].named++}if(r[st]>thr[st].q3)hit[st].q3++;if(r[st]>thr[st].q99)hit[st].q99++}}}
  res[a]={tot,...Object.fromEntries(stats.map(st=>[st,[hit[st].max/tot,hit[st].q3/tot,hit[st].q99/tot,hit[st].named/Math.max(1,hit[st].max)].map(x=>Math.round(100*x))]))}}
out.res=res;out.thr=Object.fromEntries(stats.map(s=>[s,+thr[s].max.toFixed(2)]));
console.log(key.padEnd(18),'n',String(n).padStart(5),'scored',String(obs.length).padStart(5),'rare/person',out.rare.toFixed(0).padStart(4),'| a: best/diff/z/gap  %>max (>4th, >q99, named)');
for(const a of Object.keys(res).sort((x,y)=>y-x))console.log('   a='+a,'n='+res[a].tot,stats.map(st=>`${st} ${res[a][st][0]} (${res[a][st][1]},${res[a][st][2]},${res[a][st][3]})`).join(' | '));
import fs from 'fs';fs.mkdirSync('mix1',{recursive:true});fs.writeFileSync(`mix1/${key.replace(':','_')}.json`,JSON.stringify(out));
