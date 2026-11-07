// Ensemble study: do extra algorithms catch more hidden errors at the same false-flag rate?
import fs from 'fs';import * as L from './load2.mjs';
import {analyse,distances,braycurtis,scorable,rng} from '../../swap/core.mjs';
const OUT='/tmp/claude-0/-home-claude/d674b933-5c36-5fd7-9c63-4fe2a591bf95/scratchpad/atlas/ens/';fs.mkdirSync(OUT,{recursive:true});
const all={...L.VALID,...L.NEW,...L.BOREN};for(const n of L.ploegNiches())all['ploeg_'+n.replace(/\W+/g,'_')]=()=>L.ploeg(n);
const key=process.argv[2];let ds=all[key]();
if(!/relative/i.test(ds.type||'')){const keep=[];ds.data.samples.forEach((s,i)=>{if(s.val.reduce((a,b)=>a+b,0)>=1000)keep.push(i)});
  ds={...ds,data:{feats:ds.data.feats,samples:keep.map(i=>ds.data.samples[i])},lab:keep.map(i=>ds.lab[i]),time:keep.map(i=>ds.time[i])}}
const data=ds.data,lab0=ds.lab,n=data.samples.length,S=data.samples,E=1e-6;
const M=distances(data),B=braycurtis(S);
const time=(ds.time||[]).map(v=>{const x=typeof v==='number'?v:Date.parse(v);return Number.isFinite(x)?x:NaN});const hasTime=time.length===n&&time.filter(Number.isFinite).length>n*.8;
// neighbours by fingerprint distance, once
const K=30,NB=[];for(let i=0;i<n;i++){const o=[];for(let j=0;j<n;j++)if(j!==i)o.push([M.J[i*n+j],j]);o.sort((a,b)=>a[0]-b[0]);NB.push(o.slice(0,K).map(x=>x[1]))}
// rare features: solid presence (2+ reads) in <=5% of samples (and >=2 samples)
const prev=new Map();S.forEach(s=>{for(let k=0;k<s.idx.length;k++)if(s.val[k]>=2)prev.set(s.idx[k],(prev.get(s.idx[k])||0)+1)});
const rareMax=Math.max(2,Math.floor(.05*n));const RS=S.map(s=>{const r=[];for(let k=0;k<s.idx.length;k++)if(s.val[k]>=2&&prev.get(s.idx[k])>=2&&prev.get(s.idx[k])<=rareMax)r.push(s.idx[k]);return r});
const holders=new Map();RS.forEach((r,i)=>{for(const f of r){let h=holders.get(f);if(!h){h=[];holders.set(f,h)}h.push(i)}});

function singleScore(D,lab){const s=new Float64Array(n);for(let i=0;i<n;i++){let own=Infinity,oth=Infinity;for(let j=0;j<n;j++){if(j===i)continue;const d=D[i*n+j];if(lab[j]===lab[i]){if(d<own)own=d}else if(d<oth)oth=d}s[i]=Number.isFinite(own)?Math.log((own+E)/(oth+E)):NaN}return s}
function knnScore(lab){const s=new Float64Array(n);for(let i=0;i<n;i++){const v=new Map();let own=0;NB[i].slice(0,7).forEach((j,r)=>{const w=1/(r+1);if(lab[j]===lab[i])own+=w;else v.set(lab[j],(v.get(lab[j])||0)+w)});let best=0;for(const x of v.values())if(x>best)best=x;s[i]=Math.log((best+.05)/(own+.05))}return s}
function rareScore(lab){const s=new Float64Array(n);for(let i=0;i<n;i++){const c=new Map();let own=0;for(const f of RS[i]){const seen=new Set();for(const j of holders.get(f)){if(j===i)continue;const l=lab[j];if(seen.has(l))continue;seen.add(l);if(l===lab[i])own++;else c.set(l,(c.get(l)||0)+1)}}let best=0;for(const x of c.values())if(x>best)best=x;s[i]=Math.log((best+1)/(own+1))}return s}
function timeScore(lab){const s=new Float64Array(n).fill(0);if(!hasTime)return s;for(let i=0;i<n;i++){let oth=Infinity,cand=[];for(let j=0;j<n;j++){if(j===i)continue;const d=M.J[i*n+j];if(lab[j]===lab[i]){if(Number.isFinite(time[j]))cand.push([Math.abs(time[j]-time[i]),d])}else if(d<oth)oth=d}
  if(!cand.length||!Number.isFinite(time[i])){s[i]=0;continue}cand.sort((a,b)=>a[0]-b[0]);const adj=cand.slice(0,2);const dA=adj.reduce((a,x)=>a+x[1],0)/adj.length;s[i]=Math.log((dA+E)/(oth+E))}return s}

const METHODS=['swapseq','fpJ','aitA','bray','knn','rare','time'];
function features(lab,seed){const r=analyse(data,lab,{seed},M);
  return {swapseq:r.rows.map(x=>Math.max(x.prob||0,.999*(x.broad||0))),strictFlag:r.rows.map(x=>x.flag),checkLvl:r.rows.map(x=>x.level!=='ok'),
    pairRaw:r.rows.map(x=>x.pair),fpJ:singleScore(M.J,lab),aitA:singleScore(M.A,lab),bray:singleScore(B,lab),knn:knnScore(lab),rare:rareScore(lab),time:timeScore(lab)}}
// hidden errors as in the Atlas
const ok=scorable(lab0),cnt=new Map();lab0.forEach(l=>cnt.set(l,(cnt.get(l)||0)+1));
const idx=[...Array(n).keys()].filter(i=>ok[i]&&cnt.get(lab0[i])>=3);
const small=n<300,NS=small?1:Math.max(3,Math.round(n*.003)),NO=small?1:Math.max(2,Math.round(n*.002));
if(idx.length<4*(NS+NO)){fs.writeFileSync(OUT+key+'.json',JSON.stringify({key,skipped:true}));console.log(key,'skipped');process.exit(0)}
const REPS=+(process.env.REPS||10),reps=[];
for(let rep=0;rep<REPS;rep++){const R=rng(77000+rep*13),lab=lab0.slice(),used=new Set(),pos=new Set();
  const pick=()=>{for(;;){const i=idx[Math.floor(R()*idx.length)];if(!used.has(i))return i}};let sw=0,on=0;
  while(sw<NS){const a=pick(),b=pick();if(a===b||lab[a]===lab[b])continue;[lab[a],lab[b]]=[lab[b],lab[a]];used.add(a);used.add(b);pos.add(a);pos.add(b);sw++}
  while(on<NO){const a=pick(),b=pick();if(a===b||lab[a]===lab[b])continue;lab[a]=lab[b];used.add(a);used.add(b);pos.add(a);on++}
  const F=features(lab,rep+3);const sc=scorable(lab);reps.push({F,pos:[...pos],sc})}
// rows: [features..., y] for scorable samples
const FE=['fpJ','aitA','bray','knn','rare','time','pairRaw'];
const rowsOf=rp=>{const X=[],y=[],sw=[],ss=[];for(let i=0;i<n;i++){if(!rp.sc[i])continue;X.push(FE.map(f=>{const v=rp.F[f][i];return Number.isFinite(v)?v:0}));y.push(rp.pos.includes(i)?1:0);sw.push(rp.F.swapseq[i]);ss.push([rp.F.strictFlag[i],rp.F.checkLvl[i]])}return {X,y,sw,ss}};
function fitLogit(X,y,lam=1){const p=X[0].length+1,N=X.length,n1=y.reduce((a,b)=>a+b,0),n0=N-n1;const mu=Array(p-1).fill(0),sd=Array(p-1).fill(0);
  X.forEach(r=>r.forEach((v,k)=>mu[k]+=v/N));X.forEach(r=>r.forEach((v,k)=>sd[k]+=(v-mu[k])**2/N));for(let k=0;k<p-1;k++)sd[k]=Math.sqrt(sd[k])||1;
  const Z=X.map(r=>[1,...r.map((v,k)=>(v-mu[k])/sd[k])]),w=y.map(v=>v?.5/n1:.5/n0);let b=Array(p).fill(0);
  for(let it=0;it<50;it++){const g=Array(p).fill(0),H=Array.from({length:p},()=>Array(p).fill(0));
    for(let i=0;i<N;i++){const z=Z[i];let e=0;for(let k=0;k<p;k++)e+=b[k]*z[k];const pr=1/(1+Math.exp(-e));const r=w[i]*(y[i]-pr),q=w[i]*pr*(1-pr);for(let k=0;k<p;k++){g[k]+=r*z[k];for(let l=0;l<p;l++)H[k][l]+=q*z[k]*z[l]}}
    for(let k=1;k<p;k++){g[k]-=lam*1e-3*b[k];H[k][k]+=lam*1e-3}for(let k=0;k<p;k++)H[k][k]+=1e-9;
    // solve H d = g
    const A=H.map((r,i)=>[...r,g[i]]);for(let c=0;c<p;c++){let m=c;for(let r=c+1;r<p;r++)if(Math.abs(A[r][c])>Math.abs(A[m][c]))m=r;[A[c],A[m]]=[A[m],A[c]];for(let r=0;r<p;r++){if(r===c)continue;const f=A[r][c]/A[c][c];for(let k=c;k<=p;k++)A[r][k]-=f*A[c][k]}}
    const d=A.map((r,i)=>r[p]/r[i]);let ch=0;for(let k=0;k<p;k++){b[k]+=d[k];ch+=Math.abs(d[k])}if(ch<1e-8)break}
  return x=>{let e=b[0];for(let k=0;k<p-1;k++)e+=b[k+1]*(x[k]-mu[k])/sd[k];return e}}
const half=Math.floor(REPS/2),tr=reps.slice(0,half).map(rowsOf),te=reps.slice(half).map(rowsOf);
const cat=a=>({X:a.flatMap(r=>r.X),y:a.flatMap(r=>r.y),sw:a.flatMap(r=>r.sw),ss:a.flatMap(r=>r.ss)});
const TR=cat(tr),TE=cat(te);
const stackAll=fitLogit(TR.X,TR.y),stackNew=fitLogit(TR.X.map((r,i)=>[...r,TR.sw[i]]),TR.y);
// catch rate at a false-flag rate, threshold chosen on the TRAINING negatives (as a user would)
function at(scTr,scTe,fpr){const neg=scTr.filter((v,i)=>!TR.y[i]).sort((a,b)=>b-a);const k=Math.floor(fpr*neg.length);const thr=k<1?neg[0]+1e-12:neg[k-1];
  let tp=0,P=0,fp=0,N=0;scTe.forEach((v,i)=>{if(TE.y[i]){P++;if(v>thr)tp++}else{N++;if(v>thr)fp++}});return {catch:tp/P,falsePer1000:1000*fp/N}}
const col=(D,k)=>D.X.map(r=>r[k]);
const res={key,name:ds.name,n,hasTime,positivesTest:TE.y.reduce((a,b)=>a+b,0),
  swapseqActual:{flag:TE.ss.filter((s,i)=>TE.y[i]&&s[0]).length/TE.y.filter(Boolean).length,check:TE.ss.filter((s,i)=>TE.y[i]&&s[1]).length/TE.y.filter(Boolean).length,
    falseFlagPer1000:1000*TE.ss.filter((s,i)=>!TE.y[i]&&s[0]).length/TE.y.filter(v=>!v).length,falseCheckPer1000:1000*TE.ss.filter((s,i)=>!TE.y[i]&&s[1]).length/TE.y.filter(v=>!v).length},methods:{}};
for(const fpr of [.0008,.002,.007,.01]){const o={};
  o.swapseq=at(TR.sw,TE.sw,fpr);FE.forEach((f,k)=>{o[f]=at(col(TR,k),col(TE,k),fpr)});
  o.stackAll=at(TR.X.map(stackAll),TE.X.map(stackAll),fpr);
  o.stackPlusSwapseq=at(TR.X.map((r,i)=>stackNew([...r,TR.sw[i]])),TE.X.map((r,i)=>stackNew([...r,TE.sw[i]])),fpr);
  res.methods['fpr'+fpr]=o}
fs.writeFileSync(OUT+key+'.json',JSON.stringify(res,null,1));
const f=v=>String(Math.round(v*100)).padStart(3);const m=res.methods['fpr0.002'],m2=res.methods['fpr0.01'];
console.log(key.padEnd(30),'n',n,'| @0.2%: sw',f(m.swapseq.catch),'J',f(m.fpJ.catch),'A',f(m.aitA.catch),'BC',f(m.bray.catch),'kNN',f(m.knn.catch),'rare',f(m.rare.catch),'time',f(m.time.catch),'STACK',f(m.stackAll.catch),'STACK+sw',f(m.stackPlusSwapseq.catch),'| @1%: sw',f(m2.swapseq.catch),'STACK+sw',f(m2.stackPlusSwapseq.catch));
