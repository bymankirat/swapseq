// Second opinion: the rank rule of Find-abnormality (Ye et al. 2026, MIT licence), written again on a log Bray-Curtis
// distance. Question: when SwapSeq's check list and the rank rule agree, is that strong enough to flag?
import * as C from '../../method/core.mjs';import {get} from '../robustness/loadany.mjs';import {load} from '../ye2026_comparison/loadye.mjs';import fs from 'fs';
const key=process.argv[2],REPS=+(process.env.REPS||10);
const ds=key.startsWith('ye:')?load(key.slice(3)):get(key),data=ds.data,lab=ds.lab,n=lab.length,M=C.distances(data),S=data.samples;
function logBC(S){const n=S.length,P=S.map(s=>{let t=0;for(const v of s.val)t+=v;const V=Float64Array.from(s.val,v=>Math.log(v/t*1e6+1));let sum=0;for(const v of V)sum+=v;return {I:s.idx,V,sum}}),D=new Float32Array(n*n);
  for(let a=0;a<n;a++){const A=P[a];for(let b=a+1;b<n;b++){const B=P[b];let i=0,j=0,mn=0;while(i<A.I.length&&j<B.I.length){if(A.I[i]===B.I[j]){mn+=Math.min(A.V[i],B.V[j]);i++;j++}else if(A.I[i]<B.I[j])i++;else j++}D[a*n+b]=D[b*n+a]=A.sum+B.sum>0?1-2*mn/(A.sum+B.sum):1}}return D}
const D=process.env.RD==="logbc"?logBC(S):process.env.RD==="A"?M.A:M.J; // the tool uses its own soft-Jaccard distance
// rank of every sample in every other sample's neighbour list (0 = nearest)
const rank=new Int32Array(n*n);for(let a=0;a<n;a++){const o=[...Array(n).keys()].filter(j=>j!==a).sort((x,y)=>D[a*n+x]-D[a*n+y]);o.forEach((j,r)=>rank[a*n+j]=r)}
// a person is in doubt when some pair of their samples are far in both directions (rank beyond 5% of samples);
// within such a person, samples not joined to the largest group by a close link (rank < 10 either way) are named
function rankRule(L){const by=new Map();L.forEach((l,i)=>{if(!by.has(l))by.set(l,[]);by.get(l).push(i)});const out=new Set(),lim=Math.floor(.05*n);
  for(const [,ix] of by){if(ix.length<2)continue;let bad=false;for(let a=0;a<ix.length&&!bad;a++)for(let b=a+1;b<ix.length;b++)if(rank[ix[a]*n+ix[b]]>lim&&rank[ix[b]*n+ix[a]]>lim){bad=true;break}
    if(!bad)continue;const par=new Map(ix.map(i=>[i,i])),find=x=>{while(par.get(x)!==x)x=par.get(x);return x};
    for(let a=0;a<ix.length;a++)for(let b=a+1;b<ix.length;b++)if(rank[ix[a]*n+ix[b]]<10||rank[ix[b]*n+ix[a]]<10)par.set(find(ix[a]),find(ix[b]));
    const size=new Map();for(const i of ix){const r=find(i);size.set(r,(size.get(r)||0)+1)}let big=null;for(const [r,c] of size)if(big==null||c>size.get(big))big=r;
    for(const i of ix)if(find(i)!==big)out.add(i)}return out}
const cnt=new Map();lab.forEach(l=>cnt.set(l,(cnt.get(l)||0)+1));const pool=[...Array(n).keys()].filter(i=>cnt.get(lab[i])>=3);
if(pool.length<20){console.log(key,'too few');process.exit(0)}
const NP=n<300?1:Math.max(3,Math.round(n*.003));const A={flag:[0,0],check:[0,0],checkAgree:[0,0],rule:[0,0],flagAgree:[0,0]};let T=0,O=0;
for(let rep=0;rep<REPS;rep++){const R=C.rng(5000+rep),L=lab.slice(),sw=new Set();let g=0;
  while(sw.size<2*NP&&g++<1e5){const a=pool[Math.floor(R()*pool.length)],b=pool[Math.floor(R()*pool.length)];if(a===b||L[a]===L[b]||sw.has(a)||sw.has(b))continue;[L[a],L[b]]=[L[b],L[a]];sw.add(a);sw.add(b)}
  T+=sw.size;O+=n-sw.size;const r=C.promotePairs(C.analyse(data,L,{seed:rep+3},M),null),rr=rankRule(L);
  const add=(k,i)=>A[k][sw.has(i)?0:1]++;for(const i of rr)add('rule',i);
  for(const x of r.rows){if(x.flag){add('flag',x.i);if(rr.has(x.i))add('flagAgree',x.i)}else if(x.level==='check'){add('check',x.i);if(rr.has(x.i))add('checkAgree',x.i)}}}
const f=k=>`${Math.round(100*A[k][0]/T)}% / ${(1000*A[k][1]/O).toFixed(1)}`;
console.log(key.padEnd(18),'hidden',String(T).padStart(3),'| rank rule',f('rule'),'| flags',f('flag'),'(agree',f('flagAgree')+')','| check only',f('check'),'(agree',f('checkAgree')+')');
fs.mkdirSync('out_'+(process.env.RD||'L'),{recursive:true});fs.writeFileSync(`out_${process.env.RD||"L"}/${key.replace(':','_')}.json`,JSON.stringify({key,n,T,O,A}));
