// SPDX-License-Identifier: GPL-3.0-or-later
// SwapSeq core 1.1.1 (method v3) bundled for R and the command line. Generated from core.mjs; do not edit by hand.
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mankirat Singh
// SwapSeq core: finds microbiome samples whose community does not match the person they are labelled as.
// Pure JS, no dependencies, shared by the validation scripts (Node) and the browser tool.
//
// Two tests, both self-calibrated on the user's own data:
//  1. Single-sample test (any mislabel): is this sample closer to another person than to the rest of its own person?
//     Distance = soft Jaccard on depth-matched expected presence of ASVs/OTUs (the strongest fingerprint in our benchmarks).
//  2. Pair test (two tubes swapped): would swapping the labels of sample i and sample j make both fit better?
//     Distance = mean of Jaccard and Aitchison (CLR) evidence. Also names the probable partner.
// Calibration: label swaps are planted in the user's own table, so the tool learns what a swap looks like
// in that particular study, then turns each score into a probability and reports the study's power.

// ---------- samples ----------
// A sample is {id, idx:Int32Array (sorted feature indices), val:Float64Array (counts)}.
function fromDense(featIds, sampIds, getCount){
  const S=[];
  for(let j=0;j<sampIds.length;j++){const I=[],V=[];for(let i=0;i<featIds.length;i++){const c=getCount(i,j);if(c>0){I.push(i);V.push(c)}}S.push({id:sampIds[j],idx:Int32Array.from(I),val:Float64Array.from(V)})}
  return {feats:featIds,samples:S};
}
// CSC-style sparse (BIOM 2.1 sample/matrix): indptr over samples, indices = feature rows
function fromCSC(featIds,sampIds,indptr,indices,data){
  const S=[];
  for(let j=0;j<sampIds.length;j++){const a=indptr[j],b=indptr[j+1];const I=[],V=[];
    for(let k=a;k<b;k++)if(data[k]>0){I.push(indices[k]);V.push(data[k])}
    const o=I.map((_,q)=>q).sort((x,y)=>I[x]-I[y]);
    S.push({id:sampIds[j],idx:Int32Array.from(o.map(q=>I[q])),val:Float64Array.from(o.map(q=>V[q]))})}
  return {feats:featIds,samples:S};
}

// ---------- distances (full n*n Float32Array) ----------
function braycurtis(S){ // on relative abundance: 1 - sum(min(p,q))
  const n=S.length,D=new Float32Array(n*n),P=S.map(s=>{let t=0;for(const v of s.val)t+=v;return s.val.map(v=>v/(t||1))});
  for(let a=0;a<n;a++){const ia=S[a].idx,pa=P[a];for(let b=a+1;b<n;b++){const ib=S[b].idx,pb=P[b];let i=0,j=0,m=0;
    while(i<ia.length&&j<ib.length){const x=ia[i],y=ib[j];if(x===y){m+=pa[i]<pb[j]?pa[i]:pb[j];i++;j++}else if(x<y)i++;else j++}
    D[a*n+b]=D[b*n+a]=1-m}}
  return D;
}
function jaccard(S){ // presence/absence
  const n=S.length,D=new Float32Array(n*n),T=S.map(s=>s.idx);
  for(let a=0;a<n;a++){const ia=T[a];for(let b=a+1;b<n;b++){const ib=T[b];let i=0,j=0,m=0;
    while(i<ia.length&&j<ib.length){const x=ia[i],y=ib[j];if(x===y){m++;i++;j++}else if(x<y)i++;else j++}
    const u=ia.length+ib.length-m;D[a*n+b]=D[b*n+a]=u?1-m/u:0}}
  return D;
}
function aitchison(S,nFeat,minPrev=0.05,maxFeat=3000){ // CLR (pseudocount 0.5) on features present in >= minPrev of samples
  const n=S.length,prev=new Int32Array(nFeat);for(const s of S)for(const i of s.idx)prev[i]++;
  let keep=[];for(let i=0;i<nFeat;i++)if(prev[i]>=Math.max(2,minPrev*n))keep.push(i);
  if(keep.length>maxFeat)keep=keep.sort((x,y)=>prev[y]-prev[x]).slice(0,maxFeat);keep.sort((x,y)=>x-y);
  const pos=new Int32Array(nFeat).fill(-1);keep.forEach((f,k)=>pos[f]=k);const m=keep.length;
  const X=new Float64Array(n*m);
  for(let a=0;a<n;a++){const r=X.subarray(a*m,(a+1)*m);r.fill(0.5);const s=S[a];for(let k=0;k<s.idx.length;k++){const p=pos[s.idx[k]];if(p>=0)r[p]=s.val[k]+0.5}
    let mu=0;for(let k=0;k<m;k++){r[k]=Math.log(r[k]);mu+=r[k]}mu/=m;for(let k=0;k<m;k++)r[k]-=mu}
  const D=new Float32Array(n*n);
  for(let a=0;a<n;a++){const ra=X.subarray(a*m,(a+1)*m);for(let b=a+1;b<n;b++){const rb=X.subarray(b*m,(b+1)*m);let d=0;for(let k=0;k<m;k++){const t=ra[k]-rb[k];d+=t*t}D[a*n+b]=D[b*n+a]=Math.sqrt(d)}}
  return {D,nFeat:m};
}


// Soft Jaccard on expected presence: each sample is notionally subsampled to the same depth r
// (median depth / 8), so a taxon seen c times in a sample of depth d counts as present with
// probability 1-(1-r/d)^c. Single stray reads count for little and uneven sequencing depth no
// longer makes deep samples look unlike shallow ones. Distance = 1 - sum(min)/sum(max).
function softJaccard(S,frac=1/8){
  const n=S.length,depth=S.map(s=>{let t=0;for(const v of s.val)t+=v;return t}),med=[...depth].sort((a,b)=>a-b)[n>>1],r=Math.max(1,med*frac);
  const P=S.map((s,i)=>{const q=Math.min(1,r/(depth[i]||1)),I=[],V=[];for(let k=0;k<s.idx.length;k++){const v=q>=1?1:1-Math.pow(1-q,s.val[k]);if(v>1e-3){I.push(s.idx[k]);V.push(v)}}
    let t=0;for(const v of V)t+=v;return {I:Int32Array.from(I),V:Float64Array.from(V),t}});
  const D=new Float32Array(n*n);
  for(let a=0;a<n;a++){const A=P[a];for(let b=a+1;b<n;b++){const B=P[b];let i=0,j=0,mn=0;
    while(i<A.I.length&&j<B.I.length){const x=A.I[i],y=B.I[j];if(x===y){const u=A.V[i],w=B.V[j];mn+=u<w?u:w;i++;j++}else if(x<y)i++;else j++}
    const mx=A.t+B.t-mn;D[a*n+b]=D[b*n+a]=mx>0?1-mn/mx:1}}  // sum(max) = sum(a)+sum(b)-sum(min)
  return D;
}

// ---------- scores ----------
const E=1e-6;
// For every sample and every subject: the two nearest samples of that subject [d1, j1, d2, j2].
function nearestBySubject(D,n,lab){const T=new Array(n);
  for(let i=0;i<n;i++){const m=new Map();for(let j=0;j<n;j++){if(j===i)continue;const d=D[i*n+j],l=lab[j];let t=m.get(l);if(!t){t=[Infinity,-1,Infinity,-1];m.set(l,t)}
    if(d<t[0]){t[2]=t[0];t[3]=t[1];t[0]=d;t[1]=j}else if(d<t[2]){t[2]=d;t[3]=j}}T[i]=m}
  return T}
const near=(t,ex)=>t?(t[1]===ex?t[2]:t[0]):Infinity;
// single: log(own/best other), > 0 means the sample looks more like somebody else than like its own person
function singleScores(D,n,lab){const T=nearestBySubject(D,n,lab),s=new Float64Array(n),who=new Array(n),nb=new Int32Array(n);
  for(let i=0;i<n;i++){const own=near(T[i].get(lab[i]),-1);let best=Infinity,bl=null,bj=-1;for(const [l,t] of T[i])if(l!==lab[i]&&t[0]<best){best=t[0];bl=l;bj=t[1]}
    s[i]=isFinite(own)&&isFinite(best)?Math.log((own+E)/(best+E)):NaN;who[i]=bl;nb[i]=bj}
  return {s,who,nb,T}}
// pair: max over j of [log(own_i / i-to-subject(j)) + log(own_j / j-to-subject(i))], the gain from swapping i and j
function pairScores(Dlist,n,lab){const Ts=Dlist.map(D=>nearestBySubject(D,n,lab)),s=new Float64Array(n).fill(NaN),partner=new Int32Array(n).fill(-1);
  const own=Ts.map(T=>Array.from({length:n},(_,i)=>near(T[i].get(lab[i]),-1)));
  for(let i=0;i<n;i++){let bp=-Infinity,bj=-1;
    for(let j=0;j<n;j++){if(lab[j]===lab[i])continue;let v=0,ok=true;
      for(let q=0;q<Ts.length;q++){const T=Ts[q],oi=own[q][i],oj=own[q][j],iB=near(T[i].get(lab[j]),j),jA=near(T[j].get(lab[i]),i);
        if(!isFinite(oi)||!isFinite(oj)||!isFinite(iB)||!isFinite(jA)){ok=false;break}v+=Math.log((oi+E)/(iB+E))+Math.log((oj+E)/(jA+E))}
      if(ok){v/=Ts.length;if(v>bp){bp=v;bj=j}}}
    if(bj>=0){s[i]=bp;partner[i]=bj}}
  return {s,partner}}

function scorable(lab){const c=new Map();for(const l of lab)c.set(l,(c.get(l)||0)+1);return lab.map(l=>c.get(l)>=2)}
function rng(seed){let x=seed>>>0||1;return()=>{x^=x<<13;x>>>=0;x^=x>>17;x^=x<<5;x>>>=0;return x/4294967296}}

// ---------- self-calibration by planted errors ----------
// Plant errors in the user's own labels, rescore only the planted samples, and collect the scores.
// Single test: one-way relabels (sample of A labelled as B). Pair test: reciprocal swaps.
function planted(DJ,Dpair,n,lab,nPlant,seed){
  const R=rng(seed),ok=scorable(lab),idx=[];for(let i=0;i<n;i++)if(ok[i])idx.push(i);
  const single=[],pair=[];if(idx.length<4)return {single,pair};
  const TJ=nearestBySubject(DJ,n,lab),TP=Dpair.map(D=>D===DJ?TJ:nearestBySubject(D,n,lab));
  const ownP=TP.map(T=>Array.from({length:n},(_,i)=>near(T[i].get(lab[i]),-1)));
  for(let t=0;t<nPlant;t++){const a=idx[Math.floor(R()*idx.length)],b=idx[Math.floor(R()*idx.length)];const la=lab[a],lb=lab[b];if(la===lb)continue;
    // one-way relabel: a (truly la) is labelled lb
    {const own=near(TJ[a].get(lb),-1);let best=Infinity;for(const [l,t] of TJ[a])if(l!==lb&&t[0]<best)best=t[0];
     if(isFinite(own)&&isFinite(best))single.push(Math.log((own+E)/(best+E)))}
    // reciprocal swap a<->b, scored for a: candidate partners are b and every sample outside subjects la and lb
    {let bp=-Infinity;
     const tryJ=j=>{let v=0;for(let q=0;q<TP.length;q++){const T=TP[q];let oi,oj,iB,jA;
        oi=near(T[a].get(lb),b);                       // a's own under its new label: B without b
        if(j===b){oj=near(T[b].get(la),a);iB=near(T[a].get(la),-1);jA=near(T[b].get(lb),-1)}
        else{const lc=lab[j];oj=ownP[q][j];iB=near(T[a].get(lc),j);jA=near(T[j].get(lb),b)}
        if(!isFinite(oi)||!isFinite(oj)||!isFinite(iB)||!isFinite(jA))return;v+=Math.log((oi+E)/(iB+E))+Math.log((oj+E)/(jA+E))}
       v/=TP.length;if(v>bp)bp=v};
     tryJ(b);for(let j=0;j<n;j++){if(lab[j]===la||lab[j]===lb)continue;tryJ(j)}
     if(bp>-Infinity)pair.push(bp)}}
  return {single,pair};
}
function kdeW(x,w,grid,bw){const f=new Float64Array(grid.length);let W=0;for(const v of w)W+=v;const c=1/((W||1)*bw*Math.sqrt(2*Math.PI));
  const lo=grid[0],step=grid[1]-grid[0],span=Math.ceil(6*bw/step);
  for(let i=0;i<x.length;i++){if(w[i]<=0)continue;const g0=Math.round((x[i]-lo)/step);for(let g=Math.max(0,g0-span);g<=Math.min(grid.length-1,g0+span);g++){const z=(grid[g]-x[i])/bw;f[g]+=w[i]*Math.exp(-.5*z*z)}}
  for(let g=0;g<grid.length;g++)f[g]*=c;return f}
function sd(x){const m=x.reduce((a,b)=>a+b,0)/x.length;return Math.sqrt(x.reduce((a,b)=>a+(b-m)*(b-m),0)/Math.max(1,x.length-1))}
// Empirical Bayes. All observed scores together have density g = (1-pi)*f0 + pi*f1, where f1 is what an
// error looks like (learnt from errors planted in this same table) and pi is the error rate.
// Bayes' rule then gives P(error | score) = pi*f1(score) / g(score) without ever needing f0, so real errors
// already in the table do not distort the answer. pi is either the user's prior or estimated from the
// upper tail: pi ~ #(observed >= t) / (N * P_f1(>= t)), taking the smallest value over t in the planted
// upper half (an upper bound when correct samples also reach that far). Estimation is off by default
// (opt.estimate === true switches it on); it was tested and added false flags, so the prior is used.
// Logistic version: fit log(f1/g) = a + b*score by logistic regression of planted errors (1) against all
// observed samples (0). A straight line in the log-ratio is smooth in sparse tails and automatically
// monotone, so there are no kernel spikes. P(error | score) = pi * f1/g, capped at 1.
function calibrateLogit(obs,plant,prior=0.01){
  obs=obs.filter(Number.isFinite);plant=plant.filter(Number.isFinite);if(obs.length<5||plant.length<5)return Object.assign(()=>NaN,{pi:NaN});
  const x=obs.concat(plant),y=obs.map(()=>0).concat(plant.map(()=>1)),N0=obs.length,N1=plant.length;
  const wt=y.map(v=>v?.5/N1:.5/N0); // balanced classes; the class ratio is added back below
  const mu=x.reduce((a,b)=>a+b,0)/x.length,sdv=Math.sqrt(x.reduce((a,b)=>a+(b-mu)*(b-mu),0)/x.length)||1,z=x.map(v=>(v-mu)/sdv);
  let a=0,b=0;
  for(let it=0;it<100;it++){let g0=0,g1=0,h00=0,h01=0,h11=0;
    for(let i=0;i<z.length;i++){const p=1/(1+Math.exp(-(a+b*z[i]))),r=wt[i]*(y[i]-p),w=wt[i]*p*(1-p);g0+=r;g1+=r*z[i];h00+=w;h01+=w*z[i];h11+=w*z[i]*z[i]}
    h00+=1e-9;h11+=1e-6; // tiny ridge keeps perfectly separated data finite
    const det=h00*h11-h01*h01,da=(h11*g0-h01*g1)/det,db=(h00*g1-h01*g0)/det;a+=da;b+=db;if(Math.abs(da)+Math.abs(db)<1e-9)break;if(Math.abs(b)>60){b=Math.sign(b)*60;break}}
  // with balanced weights, logit = log(f1/g)
  const P=s=>{if(!Number.isFinite(s))return NaN;if(s<=0)return 0;const lr=a+b*(s-mu)/sdv;return Math.min(1,prior*Math.exp(Math.min(50,lr)))};
  P.pi=prior;P.coef=[a,b,mu,sdv];return P;
}
// Tail-area version (no kernel smoothing, so no artefacts in sparse tails):
// P(error | score >= s) = pi * F1(>= s) / G(>= s), counted directly. F1 is read at min(s, its 90th percentile)
// so a score beyond every planted error is not penalised for being extreme. Made monotone in s.
function calibrateTail(obs,plant,prior=0.01){
  obs=obs.filter(Number.isFinite);plant=plant.filter(Number.isFinite);if(obs.length<5||plant.length<5)return Object.assign(()=>NaN,{pi:NaN});
  const O=Float64Array.from(obs).sort(),Pl=Float64Array.from(plant).sort(),N=O.length,Np=Pl.length,cap=Pl[Math.floor(.9*(Np-1))];
  const geq=(A,v)=>{let lo=0,hi=A.length;while(lo<hi){const m=(lo+hi)>>1;if(A[m]<v)lo=m+1;else hi=m}return A.length-lo};
  const raw=v=>{if(!(v>0))return 0;const g=Math.max(1,geq(O,v))/N,f=geq(Pl,Math.min(v,cap))/Np;return Math.min(1,prior*f/g)};
  // monotone: evaluate at the observed scores and carry the running maximum upward
  const xs=Array.from(O).filter(v=>v>0),ys=[];let m=0;for(const v of xs){m=Math.max(m,raw(v));ys.push(m)}
  const P=s=>{if(!Number.isFinite(s))return NaN;if(s<=0)return 0;let lo=0,hi=xs.length;while(lo<hi){const k=(lo+hi)>>1;if(xs[k]<=s)lo=k+1;else hi=k}
    return Math.max(raw(s),lo?ys[lo-1]:0)};
  P.pi=prior;return P;
}
function pava(y,w){const v=[],W=[],L=[];
  for(let i=0;i<y.length;i++){v.push(y[i]);W.push(w[i]);L.push(1);
    while(v.length>1&&v[v.length-2]>v[v.length-1]){const a=v.length-2,nw=W[a]+W[a+1];v[a]=(v[a]*W[a]+v[a+1]*W[a+1])/nw;W[a]=nw;L[a]+=L[a+1];v.pop();W.pop();L.pop()}}
  const out=[];for(let k=0;k<v.length;k++)for(let j=0;j<L[k];j++)out.push(v[k]);return out}
function calibrate(obs,plant,prior=0.01,opt={}){
  obs=obs.filter(Number.isFinite);plant=plant.filter(Number.isFinite);if(obs.length<5||plant.length<5)return Object.assign(()=>NaN,{pi:NaN});
  const all=obs.concat(plant);let lo=Infinity,hi=-Infinity;for(const v of all){if(v<lo)lo=v;if(v>hi)hi=v}lo-=1;hi+=1;const G=600,grid=Array.from({length:G},(_,k)=>lo+(hi-lo)*k/(G-1)),step=grid[1]-grid[0];
  const bw0=Math.max(2*step,1.06*sd(obs)*Math.pow(obs.length,-.2)),bw1=Math.max(2*step,1.06*sd(plant)*Math.pow(plant.length,-.2));
  const g=kdeW(obs,obs.map(()=>1),grid,bw0),f1=kdeW(plant,plant.map(()=>1),grid,bw1);
  let pi=prior;
  if(opt.estimate===true){const ps=[...plant].sort((x,y)=>x-y),N=obs.length;let best=Infinity;
    for(let q=.5;q<=.9001;q+=.05){const t=ps[Math.floor(q*(ps.length-1))];if(t<=0)continue;const F1=ps.filter(v=>v>=t).length/ps.length,G1=obs.filter(v=>v>=t).length/N;best=Math.min(best,G1/F1)}
    if(isFinite(best))pi=Math.min(0.25,Math.max(0.001,best))}
  const m0=Math.max(...g)*1e-9;
  const raw=grid.map((_,k)=>grid[k]<=0?0:Math.min(1,pi*f1[k]/(g[k]+m0)));
  // a higher score never gets a lower probability (running maximum), but only grid points where at least half a
  // sample's worth of observed data sits can raise it: in empty gaps the kernel ratio is noise, not evidence
  const sup=g.map(v=>v*obs.length*bw0*Math.sqrt(2*Math.PI)),post=new Array(G);let run=0;
  for(let k=0;k<G;k++){if(sup[k]>=.5&&raw[k]>run)run=raw[k];post[k]=sup[k]>=.5?Math.max(run,raw[k]):run}
  const P=s=>{if(!Number.isFinite(s))return NaN;if(s<=0)return 0;if(s<=lo)return post[0];if(s>=hi)return post[G-1];const t=(s-lo)/step,k=Math.floor(t);return post[k]+(post[Math.min(G-1,k+1)]-post[k])*(t-k)};
  P.pi=pi;return P;
}

// ---------- the full analysis ----------
// data: {feats, samples}; lab: subject per sample. Returns one row per sample plus study-level numbers.
// Relative abundances (fractions or percentages) are turned into pseudo-counts, with the smallest abundance in the
// table as one read: the depth-matched fingerprint and the CLR pseudocount are both built for counts, and on raw
// fractions the pseudocount of 0.5 would swamp every value.
function asCounts(S){let mn=Infinity,frac=false,maxTot=0;
  for(const s of S){let t=0;for(const v of s.val){t+=v;if(v>0&&v<mn)mn=v;if(!frac&&v!==Math.floor(v))frac=true}if(t>maxTot)maxTot=t}
  if(!frac||maxTot>100.5||!(mn>0))return {S,scale:1};
  const k=Math.min(1/mn,1e7/maxTot);
  return {S:S.map(s=>({id:s.id,idx:s.idx,val:Float64Array.from(s.val,v=>v>0?Math.max(1,Math.round(v*k)):0)})),scale:k}}
function distances(data,progress=()=>{},baseline=false){
  const ac=asCounts(data.samples),S=ac.S;if(baseline){progress('Bray-Curtis distances');const B=braycurtis(S);return {J:B,A:B,aitFeat:0}} // benchmark baseline only
progress('Fingerprint distances');const J=softJaccard(S);progress('Aitchison distances');const A=aitchison(S,data.feats.length);
  return {J,A:A.D,aitFeat:A.nFeat,empty:S.map(s=>s.idx.length===0),rescaled:ac.scale!==1};
}

// Cross-visit mode: when a visit/timepoint column is given, a sample is never compared with other
// people's samples from the same visit. Samples taken (and usually processed) together can look
// alike for reasons that have nothing to do with identity: newborns at day 1, one sequencing run,
// one extraction batch. A person's own samples come from other visits anyway, so a real mix-up
// still shows: the sample matches the true owner's samples from other visits.
function maskVisits(M,visit,lab=null){if(!visit||M.masked)return M;const n=visit.length,BIG=1e6;
  for(const D of new Set([M.J,M.A]))for(let i=0;i<n;i++){const vi=visit[i];if(vi==null||vi==='')continue;for(let j=0;j<n;j++)if(j!==i&&visit[j]===vi&&!(lab&&lab[j]===lab[i]))D[i*n+j]=BIG}
  M.masked=true;return M}
function analyse(data,lab,opt={},M=null){
  const prior=opt.prior??0.01,thr=opt.threshold??0.5,progress=opt.progress||(()=>{});
  const S=data.samples,n=S.length;M=M||distances(data,progress);if(opt.visit)maskVisits(M,opt.visit,lab);
  progress('Scoring samples');
  const one=singleScores(M.J,n,lab),pr=pairScores([M.J,M.A],n,lab);
  for(let i=0;i<n;i++)if(S[i].idx.length===0){one.s[i]=NaN;pr.s[i]=NaN}   // a sample with no reads cannot be checked
  progress('Planting test swaps in your data');
  const pl=planted(M.J,[M.J,M.A],n,lab,opt.nPlant||Math.min(1500,Math.max(300,n)),opt.seed||7);
  // Two levels, as in clinical QC. "flag" = strict kernel calibration (near-zero false alarms in our benchmarks).
  // "check" = logistic calibration, which catches more but also lists some natural look-alikes for a records check.
  const Ps=calibrate(Array.from(one.s),pl.single,prior,opt),Pp=calibrate(Array.from(pr.s),pl.pair,prior,opt);
  const Bs=calibrateLogit(Array.from(one.s),pl.single,prior),Bp=calibrateLogit(Array.from(pr.s),pl.pair,prior);
  const rows=[];
  for(let i=0;i<n;i++){const ps=Ps(one.s[i]),pp=Pp(pr.s[i]),bs=Bs(one.s[i]),bp=Bp(pr.s[i]),j=pr.partner[i];
    rows.push({i,id:S[i].id,subject:lab[i],single:one.s[i],pair:pr.s[i],pSingle:ps,pSwap:pp,prob:Math.max(ps||0,pp||0),broad:Math.max(bs||0,bp||0),bSwap:bp||0,
      looksLike:one.who[i],nearest:one.nb[i]>=0?S[one.nb[i]].id:null,partner:j>=0?S[j].id:null,partnerSubject:j>=0?lab[j]:null,scorable:Number.isFinite(one.s[i])})}
  const uncal=!Number.isFinite(Ps.pi)&&!Number.isFinite(Pp.pi);if(uncal)for(const r of rows)r.scorable=false;  // too few samples to calibrate: nothing is judged
  for(const r of rows){r.flag=r.prob>=thr;r.check=!r.flag&&(r.broad>=thr||r.prob>=.2);r.level=r.flag?'flag':r.check?'check':'ok'}
  // a swap is reported when both samples point at each other and both carry swap evidence at their level
  const pairs=[],seen=new Set(),swapEv=r=>r.flag?r.pSwap>=thr:r.bSwap>=thr;
  for(const r of rows){if(r.level==='ok'||!swapEv(r)||r.partner==null)continue;const o=rows[pr.partner[r.i]];
    if(o&&o.level!=='ok'&&pr.partner[o.i]===r.i&&swapEv(o)){const k=[r.id,o.id].sort().join('|');if(!seen.has(k)){seen.add(k);pairs.push({a:r.id,b:o.id,subjectA:r.subject,subjectB:o.subject,level:r.flag&&o.flag?'flag':'check'})}}}
  // how well planted swaps stand out from this table's own samples: share of planted swap scores above the 95th
  // percentile of all observed pair scores. Most samples are correct even when a few percent are not, so this
  // stays stable when the table already holds several errors (the measured power does not).
  const obsP=Array.from(pr.s).filter(Number.isFinite).sort((x,y)=>x-y),plP=pl.pair.filter(Number.isFinite).sort((x,y)=>x-y);
  const separation=obsP.length&&plP.length?plP.filter(v=>v>obsP[Math.floor(.95*(obsP.length-1))]).length/plP.length:NaN;
  // the same with a cut that ignores the top max(5%, 4) observed scores, so that one or two real swaps in a small
  // table do not move it; and re-identification, the share of samples that sit closer to their own person
  const cut2=obsP.length?obsP[Math.max(0,obsP.length-1-Math.max(4,Math.ceil(.05*obsP.length)))]:NaN;
  const separation2=obsP.length&&plP.length?plP.filter(v=>v>cut2).length/plP.length:NaN;
  const fin=Array.from(one.s).filter(Number.isFinite),reid=fin.length?fin.filter(v=>v<0).length/fin.length:NaN;
  // Small tables. One swap in a 30-sample table is already 7% errors, far above the expected 1%, so the calibrated
  // probability cannot reach a strong flag however clear the swap is. In tables that separate swaps very well
  // (separation2 >= 0.9 and re-identification >= 0.75), two samples that are each other's best swap partner and
  // both score like a typical planted swap (at or above the 25th percentile of planted swap scores) are reported
  // as a swap pair by that shape alone (people with three or more samples only). Tested on subsamples of 19 public tables and on 34 full tables, see validation section 12.
  if(opt.shapePairs!==false&&!uncal&&separation2>=(opt.shapeSeparation??0.9)&&reid>=(opt.shapeReid??0.75)&&plP.length){
    const cut=plP[Math.floor(.25*(plP.length-1))],cntL=new Map();for(const l of lab)cntL.set(l,(cntL.get(l)||0)+1);
    // both labels need 3 or more samples: with two per person, swapping one pair or the other pair gives the same data
    for(const r of rows){const j=pr.partner[r.i];if(j<0||!(r.pair>=cut)||cntL.get(lab[r.i])<3||cntL.get(lab[j])<3)continue;const o=rows[j];if(pr.partner[o.i]!==r.i||!(o.pair>=cut))continue;
      const k=[r.id,o.id].sort().join('|');if(seen.has(k))continue;seen.add(k);
      for(const x of [r,o])if(x.level==='ok'){x.check=true;x.level='check'}
      pairs.push({a:r.id,b:o.id,subjectA:r.subject,subjectB:o.subject,level:r.flag&&o.flag?'flag':'check',shape:true})}}
  const inPair=new Set(pairs.flatMap(p=>[p.a,p.b]));for(const r of rows)r.kind=r.level==='ok'?'ok':inPair.has(r.id)?'swap':'mismatch';
  const frac=(x,f)=>x.length?x.filter(f).length/x.length:NaN;
  const power={single:frac(pl.single,v=>Ps(v)>=thr),swap:frac(pl.pair,v=>Pp(v)>=thr),singleWithCheck:frac(pl.single,v=>Ps(v)>=.2||Bs(v)>=thr),swapWithCheck:frac(pl.pair,v=>Pp(v)>=.2||Bp(v)>=thr)};
  return {rows,pairs,power,prior,separation,separation2,reid,uncalibrated:uncal,errorRate:{single:Ps.pi,swap:Pp.pi},threshold:thr,aitFeat:M.aitFeat,nSubjects:new Set(lab).size,nScorable:rows.filter(r=>r.scorable).length,planted:pl};
}

// ---------- measured power (v1.1) ----------
// How often would SwapSeq catch an error in THIS table? Measured, not modelled: hide a few real
// errors in a copy of the user's labels, rerun the full check on the same distances with the same
// settings, and count what is caught. Each copy gets either swaps or one-way relabels (alternating),
// about 0.2% of samples and at least one, so the copy stays close to the user's own table; they are
// drawn from people with 3+ samples and never from samples already flagged. Fixed seeds, so the
// answer is repeatable. Replaces the v1.0 estimate read off the planted-score calibration, which on
// clean tables could show near 0% while the real catch rate was over 80% (SwapSeq Atlas, Nov 2026).
function measurePower(data,lab,res,M,opt={}){
  // Near-certain flags (P >= 0.9) are set aside first: they are almost surely the table's own errors,
  // and leaving them in makes every hidden test compete with them, so power reads low in tables that
  // hold errors. Borderline flags and the check list stay: many are natural look-alikes and part of
  // the table's normal variation, and removing them would make a weak table look stronger than it is.
  // Swap pairs found by shape (small tables) are set aside too, for the same reason.
  const shp=new Set((res.pairs||[]).filter(p=>p.shape).flatMap(p=>[p.a,p.b]));
  const away=opt.setAside===false?[]:res.rows.filter(r=>(r.level==='flag'&&r.prob>=.9)||shp.has(r.id)).map(r=>r.i);
  if(away.length&&!opt._reduced){const drop=new Set(away),keep=[];for(let i=0;i<data.samples.length;i++)if(!drop.has(i))keep.push(i);
    const m=keep.length,N=data.samples.length,sub=D=>{const X=new Float32Array(m*m);for(let a=0;a<m;a++){const ra=keep[a]*N;for(let b=0;b<m;b++)X[a*m+b]=D[ra+keep[b]]}return X};
    const J=sub(M.J),A=M.A===M.J?J:sub(M.A),M2={J,A,aitFeat:M.aitFeat,masked:M.masked};
    const data2={feats:data.feats,samples:keep.map(i=>data.samples[i])},lab2=keep.map(i=>lab[i]);
    const keptRows=keep.map((i,k)=>({i:k,level:res.rows[i].level}));
    return measurePower(data2,lab2,{rows:keptRows,prior:res.prior,threshold:res.threshold,nScorable:res.nScorable},M2,{...opt,_reduced:true,visit:opt.visit?keep.map(i=>opt.visit[i]):undefined})}
  const n=data.samples.length,ok=scorable(lab),reps=opt.reps??(n<300?40:n<1500?20:6);
  const base=new Set(res.rows.filter(r=>r.level!=='ok').map(r=>r.i));
  const cnt=new Map();lab.forEach(l=>cnt.set(l,(cnt.get(l)||0)+1));
  const idx=[];for(let i=0;i<n;i++)if(ok[i]&&!base.has(i)&&cnt.get(lab[i])>=3)idx.push(i);
  const k=Math.max(1,Math.round(n*.002));
  if(idx.length<8*k||new Set(idx.map(i=>lab[i])).size<3)return null;
  let s=0,s2=0,s3=0,sN=0,o=0,o2=0,oN=0,fp=0,fpN=0;
  for(let rep=0;rep<reps;rep++){
    if(opt.progress)opt.progress(`Measuring power: hidden test ${rep+1} of ${reps}`);
    const R=rng(9001+rep*7919),L=lab.slice(),used=new Set(),swaps=[],ones=[],doSwap=rep%2===0;
    const pick=()=>{for(let t=0;t<100000;t++){const i=idx[Math.floor(R()*idx.length)];if(!used.has(i))return i}return -1};
    let guard=0;
    if(doSwap)while(swaps.length<k&&guard++<10000){const a=pick(),b=pick();if(a<0||b<0||a===b||L[a]===L[b])continue;[L[a],L[b]]=[L[b],L[a]];used.add(a);used.add(b);swaps.push([a,b])}
    else while(ones.length<k&&guard++<10000){const a=pick(),b=pick();if(a<0||b<0||a===b||L[a]===L[b])continue;L[a]=L[b];used.add(a);used.add(b);ones.push(a)}
    const r=analyse(data,L,{prior:res.prior,threshold:res.threshold,seed:rep+3,visit:opt.visit},M);
    const up=new Set();for(const p of r.pairs)if(p.level==='check'&&(p.shape||r.separation>=0.7)){up.add(p.a);up.add(p.b)}
    for(const [a,b] of swaps)for(const x of [a,b]){sN++;if(r.rows[x].flag)s++;if(r.rows[x].flag||up.has(r.rows[x].id))s3++;if(r.rows[x].level!=='ok')s2++}
    for(const a of ones){oN++;if(r.rows[a].flag)o++;if(r.rows[a].level!=='ok')o2++}
    const planted=new Set(swaps.flat().concat(ones));
    for(let i=0;i<n;i++){if(planted.has(i)||base.has(i)||!r.rows[i].scorable)continue;fpN++;if(r.rows[i].flag)fp++}}
  if(!sN)return null;
  return {single:oN?o/oN:NaN,swap:s/sN,swapPairsPromoted:s3/sN,singleWithCheck:oN?o2/oN:NaN,swapWithCheck:s2/sN,falseFlagsPer1000:fpN?1000*fp/fpN:NaN,reps,tests:{swapped:sN,relabelled:oN},measured:true};
}

// ---------- reciprocal pairs (November 2026) ----------
// A swap leaves two samples that each fit the other's person. When both sit on the check list and point at each
// other, that pattern is much rarer by chance than either sample alone, so in tables where planted swaps stand
// out well (separation >= 0.7) such pairs are promoted to strong flags (pairs found by shape in small tables, see
// analyse, carry their own stricter gate). This matters most when a table holds
// several errors at once: the calibration assumes about 1% errors, and with 3% swapped samples strong flags alone
// caught few. Tested on 29 public tables, see validation section 11. In tables where swaps do not stand out it
// added false flags, so it is not applied there. When power was measured, the reported swap rate includes it.
function promotePairs(res,power,opt={}){
  const gate=res.separation>=(opt.minSeparation??0.7);
  const byId=new Map(res.rows.map(r=>[r.id,r]));let n=0;
  for(const p of res.pairs){if(p.level!=='check'||!(p.shape||gate))continue;
    const ra=byId.get(p.a),rb=byId.get(p.b);if(!ra||!rb||ra.cohab||rb.cohab||ra.downgraded||rb.downgraded)continue;  // household look-alikes stay on the check list
    for(const r of [ra,rb])if(!r.flag){r.flag=true;r.check=false;r.level='flag';r.kind='swap';r.byPair=true;n++}
    p.level='flag';p.promoted=true}
  if(power&&Number.isFinite(power.swapPairsPromoted)){power.swapBeforePromotion=power.swap;power.swap=Math.max(power.swap,power.swapPairsPromoted)}
  res.promoted=n;return res;
}

// ---------- second opinion (November 2026) ----------
// A second, much simpler rule looks at every listed sample: the rank rule of Find-abnormality (Ye et al. 2026, MIT
// licence), written again here on the soft-Jaccard distance. A person is in doubt when two of their samples are far
// from each other in both directions (each ranks beyond 5% of all samples in the other's list of neighbours); within
// such a person, a sample is named when no close link (rank under 10 either way) joins it to the person's largest
// group of samples. It uses no calibration, so it errs differently from the main check.
// Only an annotation: r.second = true (also named by the rank rule) or false, on listed samples. It never makes or
// removes a flag, because on its own the rule names far too many samples (85 per 1,000 over 29 public tables, and
// more than one in four in infant cohorts). Where the two agree, the sample is a better candidate: agreement kept
// 94-95% of listed hidden swaps and dropped 41-46% of the other listed samples (validation/second_opinion).
function secondOpinion(M,lab,res){
  const D=M&&M.J,n=lab.length;if(!D||res.uncalibrated||n<20)return res;
  const by=new Map();lab.forEach((l,i)=>{if(!by.has(l))by.set(l,[]);by.get(l).push(i)});
  const rank=(a,b)=>{const d=D[a*n+b];let r=0;for(let j=0;j<n;j++)if(j!==a&&j!==b&&D[a*n+j]<d)r++;return r},lim=Math.floor(.05*n),done=new Map();
  const judge=l=>{if(done.has(l))return done.get(l);const ix=by.get(l)||[],named=new Set();done.set(l,named);if(ix.length<2)return named;
    const m=ix.length,R=Array.from({length:m},()=>new Int32Array(m));for(let a=0;a<m;a++)for(let b=0;b<m;b++)if(a!==b)R[a][b]=rank(ix[a],ix[b]);
    let bad=false;for(let a=0;a<m&&!bad;a++)for(let b=a+1;b<m;b++)if(R[a][b]>lim&&R[b][a]>lim){bad=true;break}
    if(!bad)return named;if(m===2){named.add(ix[0]);named.add(ix[1]);return named}
    const par=ix.map((_,k)=>k),find=x=>{while(par[x]!==x)x=par[x];return x};
    for(let a=0;a<m;a++)for(let b=a+1;b<m;b++)if(R[a][b]<10||R[b][a]<10)par[find(a)]=find(b);
    const size=new Map();for(let a=0;a<m;a++){const r=find(a);size.set(r,(size.get(r)||0)+1)}let big=-1;for(const [r,c] of size)if(big<0||c>size.get(big))big=r;
    for(let a=0;a<m;a++)if(find(a)!==big)named.add(ix[a]);return named};
  let agree=0,listed=0;for(const r of res.rows){if(r.level==='ok'||!r.scorable)continue;listed++;r.second=judge(lab[r.i]).has(r.i);if(r.second)agree++}
  res.second={listed,agree};return res;
}

// ---------- wide screen (November 2026) ----------
// The same check with a wider net: the expected error rate set to 5%. Samples that would be listed then, but are not
// listed at the chosen rate, are returned as a second list ("wide screen"). In the seven shotgun tables of Ye et al.
// (2026) strong flags plus check list caught 89% of hidden swaps at 1% and 99% at 5%; on the published labels of 23
// other public tables the 5% setting listed 31 samples per 1,000 against 7.5. Nothing in the main result changes.
function wideScreen(data,lab,res,M,extra={},opt={}){
  const prior=opt.prior??.05;res.wide={prior,used:false,rows:[]};
  if(res.uncalibrated||!((res.prior??.01)<prior))return res;
  const w=promotePairs(annotate(data,lab,analyse(data,lab,{prior,seed:opt.seed,visit:opt.visit},M),extra),null),base=new Map(res.rows.map(r=>[r.id,r]));
  for(const x of w.rows){if(x.level==='ok')continue;const b=base.get(x.id);if(!b||b.level!=='ok'||!b.scorable)continue;b.wide=true;
    res.wide.rows.push({id:x.id,subject:x.subject,kind:x.kind,strong:x.level==='flag',looksLike:x.kind==='swap'?x.partnerSubject:x.looksLike,partner:x.kind==='swap'?x.partner:null,nearest:x.nearest,rare:x.rare||null,well:x.well||null,cohab:x.cohab||null})}
  res.wide.used=true;return res;
}

// ---------- duplicate check (November 2026) ----------
// Duplicate check: two samples labelled as different people whose communities are near-identical.
// That happens when one tube is sequenced twice under two names, a DNA extract is reused, or a file is
// uploaded twice; biological samples from two people are almost never that close. Works with one sample per
// person (cross-sectional studies) as well as with repeat samples.
//
// For every sample, its nearest neighbour among OTHER people's samples (soft Jaccard and Aitchison). A pair is a
// candidate when the two are each other's nearest neighbour and, on both distances, closer than almost any other
// pair of samples from different people: log distance below median - k * MAD (k = 5) of all samples'
// nearest-other-person log distances. With repeat samples, the pair must also be closer than 0.7 times the distance
// from either sample to its own person's nearest other sample. Tested on the 15 shotgun tables of Ye et al. (2026)
// with planted technical duplicates (abundances jittered, rare taxa dropped): 97% found in cross-sectional tables,
// 63% in longitudinal ones, and no pair reported that was not planted. In repeat-sample studies a sample put under
// the wrong person can also show up as a near-identical pair (7 of 105 planted mislabels), which is why the
// message names both causes.
function findDuplicates(M,lab,opt={},ids=null){
  const n=lab.length,k=opt.k??5,Ds=M.A&&M.A!==M.J?[M.J,M.A]:[M.J];
  const med=a=>{const s=Float64Array.from(a).sort();return s.length?s[s.length>>1]:NaN};
  const lim=[],nn=[],own=[],mode=opt.same??'pair',f=opt.f??0.7,empty=M.empty||null;
  for(const D of Ds){const best=new Float64Array(n).fill(Infinity),who=new Int32Array(n).fill(-1);
    for(let i=0;i<n;i++){if(empty&&empty[i])continue;for(let j=0;j<n;j++){if(j===i||lab[j]===lab[i]||(empty&&empty[j]))continue;const d=D[i*n+j];if(d<best[i]){best[i]=d;who[i]=j}}}
    const L=Array.from(best).filter(v=>isFinite(v)&&v>0).map(Math.log),m=med(L),mad=1.4826*med(L.map(v=>Math.abs(v-m)));
    const usable=L.length>=8&&mad>0;   // with fewer samples, or no spread, only exact copies (distance 0) are reported
    let same=Infinity;
    const sp=[],ow=new Float64Array(n).fill(Infinity);for(let i=0;i<n;i++){let b=Infinity;for(let j=0;j<n;j++){if(j!==i&&lab[j]===lab[i]){const d=D[i*n+j];if(d<b)b=d}}ow[i]=b;if(isFinite(b))sp.push(b)}
    if(mode==='global'&&sp.length>=10){const s=Float64Array.from(sp).sort();same=s[Math.floor(.05*(s.length-1))]}
    lim.push(usable?Math.min(Math.exp(m-k*mad),same):0);nn.push(who);own.push(ow)}
  const pairs=[],seen=new Set();
  for(let i=0;i<n;i++){const j=nn[0][i];if(j<0||nn[0][j]!==i)continue;            // mutual nearest neighbours on the fingerprint distance
    if(!Ds.every((D,q)=>D[i*n+j]<=lim[q])||!(Ds[0][i*n+j]<1))continue;
    if(mode==='pair'&&!Ds.every((D,q)=>D[i*n+j]<f*Math.min(own[q][i],own[q][j])))continue; // closer than either sample is to its own person's other samples
    const key=i<j?i+'|'+j:j+'|'+i;if(seen.has(key))continue;seen.add(key);
    const a=Math.min(i,j),b=Math.max(i,j);pairs.push({a,b,idA:ids?ids[a]:a,idB:ids?ids[b]:b,subjectA:lab[a],subjectB:lab[b],distance:Ds.map(D=>D[i*n+j]),limit:lim.slice()})}
  return {pairs,limits:lim};
}

// ---------- explaining flags ----------
// Plate wells like A1, A01, h12 -> {r, c, name}
function parseWell(w){const m=String(w??'').trim().toUpperCase().match(/^([A-P])\s*0?(\d{1,2})$/);if(!m)return null;const c=+m[2];if(c<1||c>48)return null;return {r:m[1].charCodeAt(0)-65,c:c-1,name:m[1]+String(c).padStart(2,'0')}}
// Adds evidence a person can check to every non-ok row:
//  rare   - rare ASVs (solid presence = 2+ reads, carried by <= 5% of samples) the sample shares with the
//           person it looks like vs with its own labelled person
//  cohab  - the look-alike lives in the same household (flag is then downgraded to "check")
//  wells  - plate positions of a swap pair, and whether they are neighbours or mirror positions in a column
function annotate(data,lab,res,extra={}){
  const S=data.samples,n=S.length,P=S.map(s=>{const o=[];for(let k=0;k<s.idx.length;k++)if(s.val[k]>=2)o.push(s.idx[k]);return o});
  const prev=new Map();for(const p of P)for(const f of p)prev.set(f,(prev.get(f)||0)+1);
  const rareMax=Math.max(2,Math.floor(.05*n)),byPerson=new Map(),idxOf=new Map(S.map((s,i)=>[s.id,i]));
  lab.forEach((l,i)=>{if(!byPerson.has(l))byPerson.set(l,[]);byPerson.get(l).push(i)});
  const featsOf=(person,ex)=>{const set=new Set();for(const j of byPerson.get(person)||[])if(!ex.has(j))for(const f of P[j])set.add(f);return set};
  if(extra.well&&extra.plateBoost!==false)plateBoost(res,extra.well,extra.plate,{k:extra.plateK});
  const house=extra.house,well=extra.well?extra.well.map(parseWell):null,plate=extra.plate;
  for(const r of res.rows){if(r.level==='ok')continue;const i=r.i,other=r.kind==='swap'?r.partnerSubject:r.looksLike;if(other==null)continue;
    const ex=new Set([i]);if(r.kind==='swap'&&idxOf.has(r.partner))ex.add(idxOf.get(r.partner));
    const rare=P[i].filter(f=>prev.get(f)<=rareMax),A=featsOf(r.subject,ex),B=featsOf(other,ex);
    r.rare={total:rare.length,withOwn:rare.filter(f=>A.has(f)).length,withOther:rare.filter(f=>B.has(f)).length,other,maxSamples:rareMax};
    if(house){const h=house[i];if(h){const hs=new Set((byPerson.get(other)||[]).map(j=>house[j]).filter(Boolean));if(hs.has(h)){r.cohab=h;if(r.level==='flag'){r.level='check';r.flag=false;r.check=true;r.downgraded=true}}}}
    if(well&&well[i])r.well=(plate&&plate[i]?plate[i]+' ':'')+well[i].name}
  res.shifts=findShifts(lab,res.rows,extra.well,plate);const inShift=new Set(res.shifts.flatMap(x=>x.ids));for(const r of res.rows)if(inShift.has(r.id)){r.shift=res.shifts.find(x=>x.ids.includes(r.id));if(r.level==='ok'){r.level='check';r.kind='mismatch'}}
  for(const p of res.pairs){const a=idxOf.get(p.a),b=idxOf.get(p.b);const ra=res.rows[a],rb=res.rows[b];
    if(ra.downgraded||rb.downgraded||ra.level==='check'||rb.level==='check')p.level='check';
    if(well&&well[a]&&well[b]&&(!plate||plate[a]===plate[b])){const wa=well[a],wb=well[b],d=Math.abs(wa.r-wb.r)+Math.abs(wa.c-wb.c);
      p.wells={a:wa.name,b:wb.name,plate:plate?plate[a]:null,relation:d===1?'adjacent':(wa.c===wb.c&&wa.r+wb.r===7)?'mirror':(d===2&&wa.r!==wb.r&&wa.c!==wb.c)?'diagonal':'apart'}}}
  return res;
}

// ---------- plate-aware prior (November 2026) ----------
// A swap between two tubes that sit next to each other on the plate is far more likely than a swap between two
// tubes at random positions: with n samples, a given sample has up to 4 direct neighbours but n-1 possible partners.
// When two samples are each other's best swap partner AND sit in directly neighbouring wells of the same plate
// (same row or column, one step apart; or next to each other in a plating order), their swap probability is recomputed with the expected error rate
// multiplied by k (default 5). Nothing is lowered for pairs that are not neighbours, so the check is never weaker
// than without a plate map. Needs unique well positions (per plate); otherwise it does nothing.
function plateBoost(res,wellRaw,plate,opt={}){
  const k=opt.k??5,thr=res.threshold??.5;res.plate={used:false,raised:0,k};
  if(!wellRaw||res.uncalibrated||!res.planted||!(k>1))return res;
  // positions are wells (A01..H12) or, failing that, a plating order (whole numbers: neighbours are one step apart)
  let well=wellRaw.map(parseWell);
  if(well.filter(Boolean).length<8){well=wellRaw.map(w=>{const t=String(w??'').trim(),v=+t;return t!==''&&Number.isInteger(v)?{r:0,c:v,name:t}:null});res.plate.order=true}
  const seen=new Set();
  for(const r of res.rows){const w=well[r.i];if(!w)continue;const key=(plate&&plate[r.i]!=null?plate[r.i]:'')+'|'+w.name;if(seen.has(key)){res.plate.reason='two samples share a position';return res}seen.add(key)}
  if(seen.size<8){res.plate.reason='too few positions';return res}
  const byId=new Map(res.rows.map(r=>[r.id,r])),Pk=calibrate(res.rows.map(r=>r.pair),res.planted.pair,Math.min(.5,(res.prior??.01)*k));
  if(!Number.isFinite(Pk.pi))return res;res.plate.used=true;
  const pairKey=new Set(res.pairs.map(p=>[p.a,p.b].sort().join('|')));
  for(const r of res.rows){if(!r.scorable||r.partner==null)continue;const o=byId.get(r.partner);if(!o||o.partner!==r.id||!o.scorable)continue;
    const wa=well[r.i],wb=well[o.i];if(!wa||!wb||(plate&&plate[r.i]!==plate[o.i]))continue;if(Math.abs(wa.r-wb.r)+Math.abs(wa.c-wb.c)!==1)continue;
    r.neighbour=o.id;const p=Pk(r.pair);if(!(p>(r.pSwap||0)))continue;
    r.pSwapPlate=p;if(p>r.prob)r.prob=p;
    if(!r.flag&&r.prob>=thr){r.flag=true;r.check=false;r.level='flag';r.byPlate=true;res.plate.raised++}}
  for(const r of res.rows){if(!r.byPlate)continue;const o=byId.get(r.partner);if(o.level==='ok')continue;r.kind='swap';o.kind='swap';
    const key=[r.id,o.id].sort().join('|');if(pairKey.has(key)){const p=res.pairs.find(p=>[p.a,p.b].sort().join('|')===key);if(r.flag&&o.flag)p.level='flag'}
    else{pairKey.add(key);res.pairs.push({a:r.id,b:o.id,subjectA:r.subject,subjectB:o.subject,level:r.flag&&o.flag?'flag':'check',plate:true})}}
  return res;
}

// ---------- sample-sheet shift scan ----------
// A sheet that slips by one row gives a run of consecutive tubes where each tube looks like the person
// written for its neighbour. Positions come from a well column (A01..H12, read down the columns and also
// along the rows) or from a numeric plating-order column. A tube counts when it sits closer to another person
// (score > 0) and that person is the one written on its neighbour; runs of 3+ such tubes are reported, and
// their tubes join the check list even when each alone was too weak to list.
function findShifts(lab,rows,well,plate){
  if(!well)return [];const out=[];
  const num=well.map(w=>{const t=String(w??'').trim();return /^\d+(\.\d+)?$/.test(t)?+t:null});
  const wp=well.map(parseWell),modes=num.filter(v=>v!=null).length>=wp.filter(Boolean).length?['order']:['column','row'];
  for(const mode of modes){
    const key=i=>mode==='order'?num[i]:wp[i]?(mode==='column'?wp[i].c*16+wp[i].r:wp[i].r*48+wp[i].c):null;
    const name=i=>mode==='order'?String(num[i]):wp[i].name;
    const groups=new Map();rows.forEach((r,i)=>{const k=key(i);if(k==null)return;const p=plate?plate[i]||'':'';if(!groups.has(p))groups.set(p,[]);groups.get(p).push(i)});
    for(const [p,ix] of groups){ix.sort((a,b)=>key(a)-key(b));
      for(const dir of [1,-1]){let run=[];
        const fits=q=>{const i=ix.at(q),j=ix.at(q+dir);if(j===undefined||q+dir<0)return false;const r=rows[i];return (r.level!=='ok'||r.single>0)&&r.looksLike!=null&&r.looksLike===lab[j]&&lab[j]!==lab[i]};
        for(let q=0;q<=ix.length;q++){if(q<ix.length&&fits(q))run.push(ix.at(q));else{if(run.length>=3){const ids=run.map(i=>rows[i].id);
          if(!out.some(o=>o.ids.join()===ids.join()))out.push({plate:p||null,from:name(run[0]),to:name(run.at(-1)),n:run.length,direction:dir,mode,ids})}run=[]}}}}}
  return out;
}

// ---------- sample map ----------
// Classical PCoA of a distance matrix: double-centre the squared distances and take the top two
// eigenvectors by power iteration. Returns coordinates and the share of variance on each axis.
function pcoa(D,n,iters=50){
  const B=new Float32Array(n*n),rm=new Float64Array(n);let tot=0;
  for(let i=0;i<n;i++){let s=0;for(let j=0;j<n;j++){const d=D[i*n+j];s+=d*d}rm[i]=s/n;tot+=rm[i]}tot/=n;
  for(let i=0;i<n;i++)for(let j=0;j<n;j++){const d=D[i*n+j];B[i*n+j]=-.5*(d*d-rm[i]-rm[j]+tot)}
  let trace=0;for(let i=0;i<n;i++)trace+=B[i*n+i];
  const vecs=[],vals=[],R=rng(3);
  for(let k=0;k<2;k++){let v=new Float64Array(n).map(()=>R()-.5),lam=0;
    for(let it=0;it<iters;it++){const w=new Float64Array(n);for(let i=0;i<n;i++){let s=0;const o=i*n;for(let j=0;j<n;j++)s+=B[o+j]*v[j];w[i]=s}
      for(const u of vecs){let d=0;for(let i=0;i<n;i++)d+=w[i]*u[i];for(let i=0;i<n;i++)w[i]-=d*u[i]}
      let nr=0;for(let i=0;i<n;i++)nr+=w[i]*w[i];nr=Math.sqrt(nr)||1;lam=nr;for(let i=0;i<n;i++)v[i]=w[i]/nr}
    vecs.push(v);vals.push(lam)}
  return {x:Float32Array.from(vecs[0],v=>v*Math.sqrt(vals[0])),y:Float32Array.from(vecs[1],v=>v*Math.sqrt(vals[1])),share:vals.map(l=>trace>0?l/trace:0)};
}

// ===================== strain-level fingerprints and host DNA (November 2026 update) =====================
// ---------- strain-level fingerprints ----------
// Input: the usual sparse table, where a feature is a strain marker, gene or SNV and any value > 0 means present.
// Only informative markers are used: present in at least 2 samples and in at most half of them (markers that
// everybody carries say nothing about who a sample came from). Two distances:
//   J: Jaccard on presence
//   A: Jaccard where each marker is weighted by -log(prevalence), so sharing a rare marker counts for more
function markerDistances(data,opt={}){
  const S=data.samples,n=S.length,F=data.feats.length,maxPrev=opt.maxPrev??0.5,minPrev=opt.minPrev??2;
  const prev=new Int32Array(F);for(const s of S)for(let k=0;k<s.idx.length;k++)if(s.val[k]>0)prev[s.idx[k]]++;
  const use=new Uint8Array(F);let nUse=0;for(let f=0;f<F;f++)if(prev[f]>=minPrev&&prev[f]<=maxPrev*n){use[f]=1;nUse++}
  const w=new Float64Array(F);for(let f=0;f<F;f++)if(use[f])w[f]=-Math.log(prev[f]/n);
  const P=S.map(s=>{const I=[];for(let k=0;k<s.idx.length;k++)if(s.val[k]>0&&use[s.idx[k]])I.push(s.idx[k]);I.sort((a,b)=>a-b);let t=0,tw=0;for(const f of I){t++;tw+=w[f]}return {I:Int32Array.from(I),t,tw}});
  const J=new Float32Array(n*n),A=new Float32Array(n*n);
  for(let a=0;a<n;a++){const X=P[a];for(let b=a+1;b<n;b++){const Y=P[b];let i=0,j=0,c=0,cw=0;
    while(i<X.I.length&&j<Y.I.length){const x=X.I[i],y=Y.I[j];if(x===y){c++;cw+=w[x];i++;j++}else if(x<y)i++;else j++}
    const u=X.t+Y.t-c,uw=X.tw+Y.tw-cw;J[a*n+b]=J[b*n+a]=u>0?1-c/u:1;A[a*n+b]=A[b*n+a]=uw>0?1-cw/uw:1}}
  return {J,A,aitFeat:nUse,mode:'markers',markers:nUse};
}

// ---------- host DNA ----------
// Shotgun reads from stool, swabs or saliva carry some of the person's own DNA. Read alleles at common SNPs
// identify a person however much their microbes change (newborns, antibiotics, illness). Host coverage is
// usually very low (one read at a few thousand sites), so two single samples rarely share sites. Instead, as
// in single-cell demultiplexing, each sample's reads are scored against every person's POOLED reads (all
// their other samples together, leave-one-out), and the person whose genotypes best explain the reads wins.
// Input G: one entry per sample (table order), a Map site -> [altReads, totalReads] (or a number in [0,1],
// read as one read). Population allele frequencies are estimated from all samples unless opt.af is given.
function hostCheck(G,lab,opt={}){
  // v2 (stress-tested): every person is judged against two levels read off the data itself, so the call
  // does not depend on how many reads a sample or a person happens to have:
  //   s(i,S) = mean over informative sites of [log-lik of i's read given S's pooled reads - log-lik given the
  //            population frequency]; informative = S has pooled reads at that site (leave-one-out for S = own).
  //   same level   = median of s over each sample's own label (most labels are right)
  //   unrelated level = median of s over all other people
  // Each s has a standard error from its own sites. A sample is called "other person" only when its own label
  // is clearly not the same person (likelihood ratio e^3 against a relative or a stranger), another person clearly
  // is (e^3 the other way), and that person explains the reads e^8 times better. When both the label and
  // another person match (identical twins, or one person entered under two IDs), it says so instead of flagging.
  const n=lab.length,e=opt.errorRate??0.01,minSites=opt.minSites??1500,A=opt.selfLR??3,minInf=opt.minInformative??30;
  const rd=v=>Array.isArray(v)?[+v[0],+v[1]]:[Number(v)>=.5?1:0,1];
  const people=[...new Set(lab)],pid=new Map(people.map((p,k)=>[p,k]));
  const glob=new Map();let pool;
  const reads=G.map(g=>{const o=[];if(g)for(const [k,v] of g){const [a,t]=rd(v);if(t>0)o.push([k,a,t])}return o});
  reads.forEach(r=>{for(const [k,a,t] of r){const g=glob.get(k)||[0,0];g[0]+=a;g[1]+=t;glob.set(k,g)}});
  const buildPools=skip=>{pool=people.map(()=>new Map());reads.forEach((r,i)=>{if(skip.has(i))return;const P=pool[pid.get(lab[i])];for(const [k,a,t] of r){const c=P.get(k)||[0,0];c[0]+=a;c[1]+=t;P.set(k,c)}})};
  const afOf=k=>{if(opt.af&&opt.af.has(k))return Math.min(.99,Math.max(.01,opt.af.get(k)));const g=glob.get(k);return g?Math.min(.95,Math.max(.05,(g[0]+1)/(g[1]+2))):.5};
  const pred=(k,A,T)=>{const f=afOf(k),pr=[(1-f)*(1-f),2*f*(1-f),f*f],q=[e,.5,1-e];let z=0,num=0;
    for(let g=0;g<3;g++){const w=pr[g]*Math.pow(q[g],A)*Math.pow(1-q[g],T-A);z+=w;num+=w*q[g]}return z>0?num/z:f};
  const ll=(a,t,q)=>a*Math.log(q)+(t-a)*Math.log(1-q);
  // per sample and person: mean and s.e. of the per-site gain over the population prior, plus the total (for the likelihood ratio)
  const S=new Array(n);let used=0;
  const score=skip=>{buildPools(skip);used=0;
  for(let i=0;i<n;i++){const r=reads[i];if(r.length<minSites){S[i]=null;continue}used++;
    const own=pid.get(lab[i]),stats=[];
    for(let p=0;p<people.length;p++){const P=pool[p];let m=0,sum=0,sq=0;
      for(const [k,a,t] of r){const c=P.get(k);if(!c)continue;let A=c[0],T=c[1];if(p===own&&!SKIP.has(i)){A-=a;T-=t}if(T<=0)continue;
        const d=ll(a,t,pred(k,A,T))-ll(a,t,afOf(k));m++;sum+=d;sq+=d*d}
      const mean=m?sum/m:0,sd=m>1?Math.sqrt(Math.max(1e-9,(sq-m*mean*mean)/(m-1))):1;
      stats.push({p,m,mean,se:m?sd/Math.sqrt(m):Infinity,total:sum})}
    S[i]={own,stats}}};
  // pass 1 with every sample pooled; pass 2 leaves out of the pools the samples pass 1 found clearly not to be
  // their labelled person, so a swapped-in sample no longer blurs the genotypes of the person it was put under
  let SKIP=new Set();score(SKIP);
  if(!used)return {ok:false,reason:'no sample has enough host reads',rows:G.map((_,i)=>({i,verdict:'unknown',sites:reads[i].length}))};
  const judge=()=>{
  const med=x=>{const y=x.filter(Number.isFinite).sort((a,b)=>a-b);return y.length?y[y.length>>1]:NaN};
  const ownS=[],othS=[];for(const x of S){if(!x)continue;for(const st of x.stats){if(st.m<minInf)continue;(st.p===x.own?ownS:othS).push(st.mean)}}
  const same=med(ownS),unrel=med(othS);
  if(!(same>unrel))return {ok:false,reason:'host reads do not separate people (too few reads shared between samples)',rows:G.map((_,i)=>({i,verdict:'unknown',sites:reads[i].length}))};
  const rows=[];
  for(let i=0;i<n;i++){const x=S[i];if(!x){rows.push({i,verdict:'unknown',sites:reads[i].length});continue}
    // position between the two levels: 0 = a stranger, 1 = the same person; relatives sit near 0.5
    const ok=st=>st&&st.m>=minInf,span=same-unrel,f=st=>(st.mean-unrel)/span,sf=st=>st.se/span;
    // three hypotheses for each person: the same person (1), a first-degree relative (0.5) or a stranger (0);
    // selfLR = log likelihood of "same person" against the better of the other two, from the position and its s.e.
    const selfLR=st=>{const x=f(st),v=2*sf(st)*sf(st);return (Math.min((x-.5)**2,x*x)-(x-1)**2)/v};
    const isSame=st=>ok(st)&&selfLR(st)>=A;      // clearly the same person (not a relative, not a stranger)
    const notSame=st=>ok(st)&&selfLR(st)<=-A;    // clearly not the same person
    // the person a sample really comes from: clearly not a stranger, and closer to "same" than to "relative"
    const vsStranger=st=>{const x=f(st);return (x*x-(x-1)**2)/(2*sf(st)*sf(st))};
    const likelyOwner=st=>ok(st)&&f(st)>=.75&&vsStranger(st)>=A;
    const o=x.stats.find(st=>st.p===x.own);
    let best=null;for(const st of x.stats){if(st.p===x.own||!ok(st))continue;if(!best||st.mean>best.mean)best=st}
    const llr=best&&o?best.total-o.total:NaN;let verdict='unclear',looks=null;
    if(isSame(o)&&best&&isSame(best)){verdict='identical DNA';looks=people[best.p]}
    else if((notSame(o)||!ok(o))&&best&&likelyOwner(best)&&llr>=(opt.cut??8)){verdict='other person';looks=people[best.p]}
    else if(isSame(o))verdict='same person';
    rows.push({i,verdict,llr,notSame:notSame(o),looksLike:looks??(best?people[best.p]:null),sites:reads[i].length,informative:o?o.m:0,ownPos:ok(o)?f(o):null,ownPosSE:ok(o)?sf(o):null,bestPos:ok(best)?f(best):null,bestPosSE:ok(best)?sf(best):null,all:opt.debug?x.stats.map(st=>({who:people[st.p],f:ok(st)?f(st):null,sf:ok(st)?sf(st):null,m:st.m})):undefined})}
  return {ok:true,rows,minSites,used,sameLevel:same,unrelatedLevel:unrel}};
  let out=judge();if(!out.ok||opt.passes===1)return out;
  SKIP=new Set(out.rows.filter(r=>r.notSame).map(r=>r.i));if(!SKIP.size)return out;
  score(SKIP);out=judge();out.setAside=SKIP.size;return out;
}

// Microbes and host DNA are independent witnesses. A host-DNA mismatch flags a sample; a host-DNA match is recorded (hostClears) but does not remove a microbial flag.
function addHost(res,host){
  if(!host||!host.ok)return res;
  for(const r of res.rows){const h=host.rows[r.i];r.host=h?{verdict:h.verdict,looksLike:h.verdict==='other person'?h.looksLike:null,llr:h.llr??null,sites:h.sites}:null;
    if(!h)continue;
    if(h.verdict==='other person'){r.level='flag';r.flag=true;r.check=false;if(r.kind==='ok')r.kind='mismatch';r.byHost=true;r.looksLike=h.looksLike||r.looksLike;r.prob=Math.max(r.prob||0,.99)}
    else if(h.verdict==='same person'&&r.level!=='ok'&&!r.byHost){r.hostClears=true}}
  return res;
}

// ---------- report (a self-contained HTML page for the lab's records) ----------
// o: {version,file,metaFile,group,date,type,prior,columns,rows,pairs,dups:[{a,sa,b,sb}],power,nScorable,nSubjects,host,uncalibrated,powerSkipped,
//     planted:[single scores of planted errors],map:{x,y} (one point per row, optional),methods,notes:[...]}
// rows are the rows of analyse/annotate. No scripts and no outside files: it opens anywhere and prints to PDF. Figures are inline SVG.
function reportHTML(o){
  const e=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const pct=v=>v==null||!Number.isFinite(v)?'not measured':Math.round(100*v)+'%';
  const rows=o.rows||[],flags=rows.filter(r=>r.level==='flag').sort((a,b)=>(b.prob||0)-(a.prob||0)),checks=rows.filter(r=>r.level==='check').sort((a,b)=>Math.max(b.broad||0,b.prob||0)-Math.max(a.broad||0,a.prob||0));
  const pw=o.power&&o.power.measured!==false?o.power:null,nPeople=new Set(rows.filter(r=>r.scorable).map(r=>r.subject)).size,notChecked=rows.filter(r=>!r.scorable&&!r.byHost).length;
  const few=nPeople<8,low=!pw||few||!(pw.swapWithCheck>=.25),dups=o.dups||[],prior=o.prior??.01,fp=(o.pairs||[]).filter(p=>p.level==='flag').length;
  const nOk=rows.length-flags.length-checks.length-notChecked,word=o.type==='markers'?'markers':'taxa';
  let cls,head,text;
  if(o.uncalibrated){cls='weak';head='Not checked.';text='Too few samples come from people with two or more samples, so no sample could be judged.'}
  else if(!pw&&o.powerSkipped&&!few){cls=flags.length?'warn':'weak';head=flags.length?`${flags.length} sample${flags.length>1?'s':''} flagged${checks.length?`, and ${checks.length} to check`:''}.`:`No strong flags${checks.length?`, ${checks.length} to check`:''}.`;text='Power was not measured in this run, so this report cannot say how many errors would be missed.'}
  else if(few||!pw){cls='weak';head='Too few people to check well.';text=(few?`Only ${nPeople} people have two or more samples.`:'Too few people have three or more samples to hide test errors in.')+' A clean result here proves little. Strong flags are still listed.'}
  else if(low){cls='weak';head='Low power in this dataset.';text=`Even with the check list, only about ${pct(pw.swapWithCheck)} of hidden swaps could be told apart from normal variation. A clean result here does not show that the labels are right.`}
  else if(flags.length){cls='warn';head=`${flags.length} sample${flags.length>1?'s':''} flagged${fp?`, including ${fp} likely swap${fp>1?'s':''}`:''}${checks.length?`, and ${checks.length} to check`:''}.`;text='Flags carry strong evidence. Compare each one with the laboratory records before changing any label.'}
  else if(checks.length){cls='weak';head=`No strong flags, ${checks.length} sample${checks.length>1?'s':''} to check.`;text='These look somewhat like another person. Some will be natural look-alikes, so compare them with the laboratory records.'}
  else{cls='good';head='No sample looks like it belongs to someone else.';text=`In this dataset SwapSeq would have caught about ${pct(pw.swapWithCheck)} of swapped tubes.`}
  const why=r=>{const a=[];
    if(r.kind==='swap'&&r.partner)a.push(`Swapping labels with ${e(r.partner)} (labelled ${e(r.partnerSubject)}) makes both samples fit their person.`);
    else if(r.nearest||r.looksLike)a.push(`Closest sample is ${r.nearest?e(r.nearest)+' from ':''}${e(r.looksLike)}, not one of ${e(r.subject)}'s own samples.`);
    if(r.rare&&r.rare.total)a.push(`Shares ${r.rare.withOther} of its ${r.rare.total} rare ${word} with ${e(r.rare.other)}'s other samples and ${r.rare.withOwn} with ${e(r.subject)}'s.`);
    if(r.host&&r.host.verdict==='other person')a.push(`Host DNA fits ${e(r.host.looksLike)} better than ${e(r.subject)}.`);
    else if(r.host&&r.host.verdict==='identical DNA')a.push(`Host DNA is identical to ${e(r.host.looksLike)} (twins, or one person under two IDs).`);
    else if(r.host&&r.host.verdict==='same person')a.push(`Host DNA fits ${e(r.subject)}, the labelled person.`);
    if(r.byPair)a.push('Strong because both samples of the pair point at each other\'s person.');
    if(r.byPlate)a.push('Strong because its swap partner sits in the next well, where a swap is more likely.');
    if(r.second===true)a.push('Second check (rank rule) agrees.');else if(r.second===false&&r.level==='check')a.push('Second check does not name it.');
    if(r.shift)a.push(`Part of a possible sample-sheet shift (positions ${e(r.shift.from)} to ${e(r.shift.to)}).`);
    if(r.cohab)a.push(`Same household (${e(r.cohab)}) as the look-alike${r.downgraded?', so it was moved to the check list':''}.`);
    if(r.well)a.push(`Well ${e(r.well)}.`);return a.join(' ')};
  const tab=(list,cap)=>list.length?`<div class="tw"><table><thead><tr><th>Sample</th><th>Labelled as</th><th>Looks like</th><th>Kind</th><th>Why</th><th class="box">Records checked</th></tr></thead><tbody>${list.slice(0,cap).map(r=>`<tr><td class="m">${e(r.id)}</td><td class="m">${e(r.subject)}</td><td class="m">${e(r.kind==='swap'?r.partnerSubject:r.looksLike)}</td><td>${r.byHost?'Host DNA':r.kind==='swap'?'Swap':'Mismatch'}</td><td>${why(r)}</td><td class="box"></td></tr>`).join('')}</tbody></table></div>${list.length>cap?`<p class="s">${list.length-cap} more are in the results table (.tsv).</p>`:''}`:'<p class="s">None.</p>';
  const set=[['Feature table',o.file],o.metaFile?['Metadata',o.metaFile]:null,o.group&&o.group!=='All samples'?['Group',o.group]:null,['Columns used',o.columns],['Data type',o.type==='markers'?'strain markers (presence)':'taxon counts'],['Expected error rate',prior*100+'%'],['Date',o.date],['Software',o.version]].filter(x=>x&&x[1]);
  // ----- figures -----
  const K={ok:'#8f99ad',check:'#c77700',flag:'#c0392b',none:'#d9dde5',strong:'#1f6f8b',soft:'#9cc7d8',own:'#2f5fd0',oth:'#1d8a55',ink:'#1a2233',mut:'#5b667a',line:'#d5dae3'};
  const fig=(n,title,svg,cap)=>`<figure><figcaption><b>Figure ${n}. ${e(title)}</b></figcaption>${svg}<p class="s">${cap}</p></figure>`;let fn=0;const figs={};
  // 1: all samples at a glance
  {const parts=[['Fit their own person',nOk,K.ok],['To check',checks.length,K.check],['Strong flag',flags.length,K.flag],['Not checked (one sample only)',notChecked,K.none]].filter(p=>p[1]>0),tot=parts.reduce((a,p)=>a+p[1],0)||1,W=800,gap=3;
    let free=W-gap*(parts.length-1);const w=parts.map(p=>Math.max(10,p[1]/tot*free)),sc=free/w.reduce((a,b)=>a+b,0);let x=0;
    const bars=parts.map((p,i)=>{const ww=w[i]*sc,r=`<rect x="${x.toFixed(1)}" y="0" width="${ww.toFixed(1)}" height="30" rx="4" fill="${p[2]}"/>`;x+=ww+gap;return r}).join('');
    figs.status=fig(++fn,'All samples at a glance',`<svg viewBox="0 0 ${W} 30" role="img" aria-label="Bar showing how many samples fit, are to check, are flagged">${bars}</svg><div class="lg">${parts.map(p=>`<span><i style="background:${p[2]}"></i>${e(p[0])}: <b>${p[1].toLocaleString('en')}</b></span>`).join('')}</div>`,
      'Each sample is in one group. Most samples should fit their own person. A strong flag has strong evidence of a wrong label. "To check" means it looks somewhat like another person.')}
  // whose rare taxa each listed sample carries
  {const L=flags.concat(low?[]:checks).filter(r=>r.rare&&r.rare.total).slice(0,12);
    if(L.length){const bar=(lab,k,tot,col)=>`<div class="eb"><span>${lab}</span><i><i style="width:${(100*k/tot).toFixed(1)}%;min-width:${k?3:0}px;background:${col}"></i></i><b>${k} of ${tot}</b></div>`;
      figs.rare=fig(++fn,`Whose rare ${word} does each listed sample carry?`,`<div class="ev">${L.map(r=>`<div><div><span class="m"><b>${e(r.id)}</b></span> <span class="s">${r.level==='flag'?'strong flag':'to check'}</span></div>${bar(`Shared with <b>${e(r.subject)}</b>, the person on its label`,r.rare.withOwn,r.rare.total,K.own)}${bar(`Shared with <b>${e(r.rare.other)}</b>, the person it looks like`,r.rare.withOther,r.rare.total,K.oth)}</div>`).join('')}</div>`,
        `Rare ${word} are found in few people, so they point to one person. A sample with a wrong label shares many with the person it looks like (green) and few with the person on its label (blue).`)}}
  // 2: power
  if(pw){const W=800,L=190,R=150,bw=W-L-R,row=(y,lab,strong,all)=>{const a=Math.max(0,Math.min(1,all||0)),s=Math.max(0,Math.min(a,strong||0));
      return `<text x="${L-10}" y="${y+16}" text-anchor="end" class="t">${lab}</text><rect x="${L}" y="${y}" width="${bw}" height="22" rx="4" fill="#eef1f5"/><rect x="${L}" y="${y}" width="${(a*bw).toFixed(1)}" height="22" rx="4" fill="${K.soft}"/><rect x="${L}" y="${y}" width="${(s*bw).toFixed(1)}" height="22" rx="4" fill="${K.strong}"/><text x="${L+bw+10}" y="${y+16}" class="t"><tspan font-weight="700">${pct(all)}</tspan> found</text>`};
    const has1=Number.isFinite(pw.singleWithCheck);
    figs.power=fig(++fn,'How many hidden errors this check finds in your data',`<svg viewBox="0 0 ${W} ${has1?92:58}" role="img" aria-label="Bars showing the share of hidden errors found">${row(4,'Two tubes swapped',pw.swap,pw.swapWithCheck)}${has1?row(38,'One sample mislabelled',pw.single,pw.singleWithCheck):''}${[0,.5,1].map(v=>`<text x="${L+v*bw}" y="${has1?86:52}" text-anchor="middle" class="a">${v*100}%</text>`).join('')}</svg><div class="lg"><span><i style="background:${K.strong}"></i>Found as a strong flag</span><span><i style="background:${K.soft}"></i>Found only on the check list</span></div>`,
      `SwapSeq swapped labels on purpose in copies of your data and counted how many it found again. Here: if 100 tubes had been swapped, about ${Math.round(100*pw.swapWithCheck)} would be found${Number.isFinite(pw.swap)?` (${Math.round(100*pw.swap)} as strong flags)`:''}. ${Math.round(100*pw.swapWithCheck)<100?'The rest would be missed.':''}`)}
  // 3: scores
  {const obs=rows.filter(r=>r.scorable&&Number.isFinite(r.single)).map(r=>r.single),pl=(o.planted||[]).filter(Number.isFinite);
    if(obs.length>=10){const all=obs.concat(pl).sort((a,b)=>a-b),lo=all[Math.floor(.005*all.length)],hi=all[Math.min(all.length-1,Math.ceil(.995*all.length))],span=(hi-lo)||1,B=36,W=800,H=230,L=14,R=14,T=14,Bm=44,pw_=W-L-R,ph=H-T-Bm;
      const bin=v=>Math.max(0,Math.min(B-1,Math.floor((v-lo)/span*B))),h1=new Array(B).fill(0),h2=new Array(B).fill(0);obs.forEach(v=>h1[bin(v)]+=1/obs.length);pl.forEach(v=>h2[bin(v)]+=1/pl.length);
      const mx=Math.max(...h1,...h2)||1,X=v=>L+Math.max(0,Math.min(1,(v-lo)/span))*pw_,Y=v=>T+ph*(1-v/mx),bwid=pw_/B;
      let g=h1.map((v,i)=>v?`<rect x="${(L+i*bwid+1).toFixed(1)}" y="${Y(v).toFixed(1)}" width="${(bwid-2).toFixed(1)}" height="${(T+ph-Y(v)).toFixed(1)}" fill="${K.ok}"/>`:'').join('');
      if(pl.length>=10){let d=`M${L},${T+ph}`;h2.forEach((v,i)=>{d+=`L${(L+i*bwid).toFixed(1)},${Y(v).toFixed(1)}L${(L+(i+1)*bwid).toFixed(1)},${Y(v).toFixed(1)}`});d+=`L${L+pw_},${T+ph}`;g+=`<path d="${d}" fill="none" stroke="${K.flag}" stroke-width="2"/>`}
      if(lo<0&&hi>0)g+=`<line x1="${X(0)}" x2="${X(0)}" y1="${T}" y2="${T+ph}" stroke="${K.mut}" stroke-dasharray="4 3"/><text x="${X(0)+5}" y="${T+10}" class="a">equally close to both</text>`;
      g+=`<line x1="${L}" x2="${L+pw_}" y1="${T+ph}" y2="${T+ph}" stroke="${K.line}"/>`;
      const mark=(r,c)=>Number.isFinite(r.single)?`<path d="M${X(r.single).toFixed(1)},${T+ph+3}l6,10h-12z" fill="${c}"><title>${e(r.id)}</title></path>`:'';
      g+=checks.map(r=>mark(r,K.check)).join('')+flags.map(r=>mark(r,K.flag)).join('');
      g+=`<text x="${L}" y="${H-6}" class="a">&#8592; looks like its own person</text><text x="${L+pw_}" y="${H-6}" text-anchor="end" class="a">looks like someone else &#8594;</text>`;
      figs.scores=fig(++fn,'The score of every sample',`<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Histogram of sample scores">${g}</svg><div class="lg"><span><i style="background:${K.ok}"></i>Your samples</span>${pl.length>=10?`<span><i class="ln" style="border-color:${K.flag}"></i>Samples with a label swapped on purpose</span>`:''}${flags.length?`<span><i class="tr" style="border-bottom-color:${K.flag}"></i>Strong flag</span>`:''}${checks.length?`<span><i class="tr" style="border-bottom-color:${K.check}"></i>To check</span>`:''}</div>`,
        `Each grey bar counts your samples with that score. Samples on the left look like their own person. ${pl.length>=10?'The red line shows where samples land when their label is wrong. The less the grey bars and the red line overlap, the better the check works in your study. ':''}Triangles under the line mark the listed samples.`)}}
  // maps: one small map per listed sample, with only the two people that matter in colour
  if(o.map&&o.map.x&&o.map.x.length===rows.length&&rows.length>=6){const Xs=o.map.x,Ys=o.map.y,W=400,H=280,P=18,mnx=Math.min(...Xs),mxx=Math.max(...Xs),mny=Math.min(...Ys),mxy=Math.max(...Ys),sx=v=>(P+(v-mnx)/((mxx-mnx)||1)*(W-2*P)).toFixed(1),sy=v=>(P+(1-(v-mny)/((mxy-mny)||1))*(H-2*P)).toFixed(1);
    const L=flags.concat(low?[]:checks).map(r=>[r,rows.indexOf(r),r.kind==='swap'?r.partnerSubject:r.looksLike]).filter(x=>x[2]!=null).slice(0,6);
    if(L.length){const base=rows.map((r,i)=>`<circle cx="${sx(Xs[i])}" cy="${sy(Ys[i])}" r="2.2" fill="${K.none}"/>`).join('');
      const one=([r,i,other])=>{let own='',oth='';rows.forEach((q,k)=>{if(k===i)return;if(q.subject===r.subject)own+=`<circle cx="${sx(Xs[k])}" cy="${sy(Ys[k])}" r="5" fill="${K.own}" stroke="#fff" stroke-width="1.2"/>`;else if(q.subject===other)oth+=`<rect x="${(+sx(Xs[k])-4.5).toFixed(1)}" y="${(+sy(Ys[k])-4.5).toFixed(1)}" width="9" height="9" fill="${K.oth}" stroke="#fff" stroke-width="1.2"/>`});
        return `<div><div class="mt"><span class="m"><b>${e(r.id)}</b></span> labelled <b style="color:${K.own}">${e(r.subject)}</b>, looks like <b style="color:${K.oth}">${e(other)}</b></div><svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Map for ${e(r.id)}"><rect x="1" y="1" width="${W-2}" height="${H-2}" rx="6" fill="none" stroke="${K.line}"/>${base}${own}${oth}<path d="M${sx(Xs[i])},${(+sy(Ys[i])-9).toFixed(1)}l9,9l-9,9l-9,-9z" fill="${r.level==='flag'?K.flag:K.check}" stroke="${K.ink}" stroke-width="1.5"/></svg></div>`};
      figs.map=fig(++fn,'Where each listed sample sits',`<div class="mm">${L.map(one).join('')}</div><div class="lg"><span><i class="dia" style="background:${K.flag}"></i>The listed sample</span><span><i class="dot" style="background:${K.own}"></i>Other samples of the person on its label</span><span><i style="background:${K.oth};border-radius:0"></i>Samples of the person it looks like</span><span><i class="dot" style="background:${K.none}"></i>Everyone else</span></div>`,
        `Samples with similar microbes sit close together. If the label were right, the diamond would sit among the blue dots. When it sits among the green squares instead, the sample looks like that other person. The map squeezes many differences into a flat picture, so it can be less clear than ${figs.rare?'Figure 2':'the table'}; use it to look, not to decide.${flags.length+checks.length>L.length?` Showing the first ${L.length} listed samples.`:''}`)}}
  return `<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>SwapSeq report: ${e(o.file||'')}</title><style>
:root{--ink:#1a2233;--mut:#5b667a;--line:#d5dae3;--bg:#fff;--soft:#f3f5f9;--warn:#c0392b;--ok:#1d6b3a;--weak:#8a5a00}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 system-ui,-apple-system,"Segoe UI",Roboto,Arial,sans-serif}
main{max-width:860px;margin:0 auto;padding:28px max(16px,env(safe-area-inset-right,0px)) 40px max(16px,env(safe-area-inset-left,0px))}
h1{font-size:24px;margin:0 0 2px}h2{font-size:17px;margin:28px 0 8px;padding-top:14px;border-top:1px solid var(--line)}p{margin:6px 0}.s{color:var(--mut);font-size:13px}
.m{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:12.5px;overflow-wrap:anywhere}
dl.set{display:grid;grid-template-columns:max-content 1fr;gap:3px 16px;margin:12px 0}dl.set dt{color:var(--mut)}dd{margin:0;overflow-wrap:anywhere}
.v{border:1px solid var(--line);border-left:5px solid var(--mut);border-radius:6px;padding:10px 14px;margin:14px 0;background:var(--soft);font-size:15px}.v.warn{border-left-color:var(--warn)}.v.good{border-left-color:var(--ok)}.v.weak{border-left-color:var(--weak)}
.n{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:14px 0}.n div{border:1px solid var(--line);border-radius:6px;padding:8px 10px}.n b{display:block;font-size:20px}.n span{color:var(--mut);font-size:12px}
figure{margin:16px 0;padding:12px 14px;border:1px solid var(--line);border-radius:6px;break-inside:avoid}figcaption{margin-bottom:10px}svg{display:block;width:100%;height:auto}
svg .t{font:13px system-ui,-apple-system,"Segoe UI",Roboto,Arial,sans-serif;fill:var(--ink)}svg .a{font:12px system-ui,-apple-system,"Segoe UI",Roboto,Arial,sans-serif;fill:var(--mut)}
.lg{display:flex;flex-wrap:wrap;gap:4px 18px;margin:8px 0 4px;font-size:13px}.lg i{display:inline-block;width:12px;height:12px;border-radius:3px;margin-right:6px;vertical-align:-1px}.lg i.dot{border-radius:50%}.lg i.dia{border-radius:0;transform:rotate(45deg) scale(.8)}
.ev,.mm{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,280px),1fr));gap:12px}.ev>div{border:1px solid var(--line);border-radius:6px;padding:9px 11px;min-width:0;break-inside:avoid}.mm>div{min-width:0;break-inside:avoid}.mt{font-size:13px;margin-bottom:5px}
.eb{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:2px 10px;align-items:center;margin-top:7px;font-size:13px}.eb>span{grid-column:1/3;color:var(--mut)}.eb>i{display:block;height:12px;border-radius:4px;background:#eef1f5;overflow:hidden}.eb>i>i{display:block;height:100%;border-radius:4px}.lg i.ln{height:0;width:16px;border-top:2px solid;border-radius:0;vertical-align:4px;background:none}.lg i.tr{width:0;height:0;border:6px solid transparent;border-top:0;border-bottom:10px solid;border-radius:0;background:none}
.tw{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{border:1px solid var(--line);padding:5px 7px;text-align:left;vertical-align:top}th{background:var(--soft);font-weight:600}.box{width:92px}
ol,ul{margin:6px 0;padding-left:22px}li{margin:3px 0}dl.w{margin:6px 0}dl.w dt{font-weight:600;margin-top:8px}
.sign{display:grid;grid-template-columns:1fr 1fr;gap:30px;margin-top:34px}.sign div{border-top:1px solid var(--ink);padding-top:4px;color:var(--mut);font-size:12px}
footer{margin-top:26px;color:var(--mut);font-size:12px;border-top:1px solid var(--line);padding-top:10px}
@media(max-width:600px){.n{grid-template-columns:1fr 1fr}.sign{grid-template-columns:1fr}}
@media print{main{max-width:none;padding:0}h2{break-after:avoid}tr{break-inside:avoid}*{-webkit-print-color-adjust:exact;print-color-adjust:exact}@page{margin:16mm}}
</style></head><body><main>
<h1>Sample identity report</h1><p class="s">SwapSeq check for swapped and mislabelled samples</p>
<dl class="set">${set.map(([k,v])=>`<dt>${e(k)}</dt><dd>${e(v)}</dd>`).join('')}</dl>
<h2>1. Result</h2>
<div class="v ${cls}"><b>${e(head)}</b> ${e(text)}</div>
<div class="n"><div><b>${(o.nScorable??0).toLocaleString('en')}</b><span>samples checked, from ${o.nSubjects??nPeople} people</span></div><div><b>${flags.length}</b><span>strong flags${fp?` (${fp} swap pair${fp>1?'s':''})`:''}</span></div><div><b>${checks.length}</b><span>on the check list</span></div><div><b>${pw?pct(pw.swapWithCheck):'n/a'}</b><span>of hidden swaps found in this dataset</span></div></div>
${(o.notes||[]).map(t=>`<p class="s">${e(t)}</p>`).join('')}
${figs.status||''}
<h2>2. Samples to look at</h2>
<p><b>Strong flags</b> (${flags.length})</p>${tab(flags,200)}
<p style="margin-top:14px"><b>Check list</b> (${checks.length})</p>${low&&checks.length?`<p class="s">Power is low here, so most of these are natural look-alikes.</p>`:''}${tab(checks,200)}
${figs.rare||''}
${o.wide&&o.wide.used?`<p style="margin-top:14px"><b>Wide screen</b> (${o.wide.rows.length})</p><p class="s">The same check with a wider net (expected error rate ${o.wide.prior*100}% instead of ${prior*100}%). It finds more real errors and lists more look-alikes, so these are weaker leads than the lists above.</p>${o.wide.rows.length?`<div class="tw"><table><thead><tr><th>Sample</th><th>Labelled as</th><th>Looks like</th><th>Kind</th><th class="box">Records checked</th></tr></thead><tbody>${o.wide.rows.slice(0,60).map(w=>`<tr><td class="m">${e(w.id)}</td><td class="m">${e(w.subject)}</td><td class="m">${e(w.looksLike)}</td><td>${w.kind==='swap'?'Swap with '+e(w.partner):'Mismatch'}</td><td class="box"></td></tr>`).join('')}</tbody></table></div>${o.wide.rows.length>60?`<p class="s">${o.wide.rows.length-60} more are in the results table (column wide_screen).</p>`:''}`:'<p class="s">None.</p>'}`:''}
${dups.length?`<p style="margin-top:14px"><b>Possible duplicates</b> (${dups.length})</p><p class="s">Samples labelled as different people with almost the same microbes: one tube or file under two names, or a sample under the wrong person.</p><ul>${dups.slice(0,50).map(d=>`<li><span class="m">${e(d.a)}</span> (labelled ${e(d.sa)}) and <span class="m">${e(d.b)}</span> (labelled ${e(d.sb)})</li>`).join('')}</ul>`:''}
${o.host?`<p style="margin-top:14px"><b>Host DNA</b></p><p>${o.host.ok?`${(o.host.used||0).toLocaleString('en')} samples had enough host reads to check; ${rows.filter(r=>r.byHost).length} flagged because their DNA fits another person.`:`Not used: ${e(o.host.reason)}.`}</p>`:''}
<h2>3. How much to trust this result</h2>
${figs.power||`<p>${o.powerSkipped?'Power was not measured in this run.':'Power could not be measured in this dataset: too few people have three or more samples.'} Without it, a clean result does not show that the labels are right.</p>`}
${figs.scores||''}${figs.map||''}
<h2>4. How SwapSeq checks a sample</h2><ol>
<li>The microbes of one person stay more alike over time than the microbes of two different people. So a person's samples work a little like a fingerprint.</li>
<li>For each sample, SwapSeq measures how close it is to the other samples with the same label, and how close it is to everyone else.</li>
<li>If a sample is much closer to another person than to its own, it is listed. If two samples would both fit after trading labels, SwapSeq calls it a swap.</li>
<li>To learn what "much closer" means in your study, SwapSeq swaps labels on purpose in copies of your data. This sets the cut-off and shows how many real errors it would find.</li></ol>
<h2>5. Words used in this report</h2><dl class="w">
<dt>Strong flag</dt><dd>Strong evidence that the label is wrong. Check the lab records.</dd>
<dt>Check list</dt><dd>Weaker evidence. Some are real errors, some are people who just look alike.</dd>
<dt>Swap</dt><dd>Two samples that each fit the other one's person. Usually two tubes or two rows traded places.</dd>
<dt>Mismatch</dt><dd>One sample that fits another person, with no partner found.</dd>
<dt>Share of hidden swaps found</dt><dd>How many swaps SwapSeq found again after hiding them in your data. If this is low, a clean result means little.</dd>
<dt>Expected error rate</dt><dd>How common label errors are assumed to be before looking (${prior*100}% here).</dd></dl>
<h2>6. What to do next</h2><ol>
<li>Find each listed sample in the lab records: plate map, sample sheet, freezer log.</li>
<li>Look for a reason: tubes next to each other, the same day, a shifted row or column, a typing error.</li>
<li>If the records confirm an error, correct the label or remove the sample. Write the decision in the last column of the table.</li>
<li>If the records show the label is right, keep the sample and note it as a natural look-alike.</li>
<li>Run SwapSeq again after any change.</li></ol>
<h2>7. Limits</h2><ul>
<li>A flag is a reason to check lab records. It is not proof, and no label should be changed on this report alone.</li>
<li>${pw?`About ${pct(1-pw.swapWithCheck)} of real swaps would be missed in this dataset.`:'The share of errors missed is not known for this dataset.'}</li>
<li>People who live together, twins and transplant pairs can look alike without any error.</li>
<li>A person with only one sample cannot be checked${notChecked?` (${notChecked} sample${notChecked>1?'s':''} here)`:''}.</li>
<li>SwapSeq mostly misses a sample that is a mixture of two people, and it finds few errors in newborns and infants.</li>
<li>Research use only. Not for diagnosis.</li></ul>
${o.methods?`<h2>8. Text for your Methods section</h2><p>${e(o.methods)}</p>`:''}
<div class="sign"><div>Records checked by</div><div>Date</div><div>Outcome (labels corrected, samples removed, or no action)</div><div>Signature</div></div>
<footer>${e(o.version)}. Singh M. (2026) SwapSeq. https://doi.org/10.5281/zenodo.23963612. The data were analysed on the user's own computer.</footer>
</main></body></html>`;
}
// ---- bridge: called from R (V8 package) and from the command-line tool ----
// input: {features:[...], samples:[{id, idx:[...0-based feature index], val:[...counts]}], subject:[...], prior,
//         house, well, plate,
//         type: "counts" (default: taxa or ASV counts) | "markers" (strain markers, genes or SNVs; any value > 0 = present),
//         host: optional, one entry per sample: {site: [altReads, totalReads], ...} from host reads at common SNPs}
function swapseq_run(json){
  const input0=typeof json==='string'?JSON.parse(json):json;
  if(!input0||!Array.isArray(input0.samples)||!Array.isArray(input0.subject))throw new Error('swapseq: need samples and one subject per sample');
  if(input0.subject.length!==input0.samples.length)throw new Error(`swapseq: ${input0.samples.length} samples but ${input0.subject.length} subject labels`);
  const prior0=input0.prior==null?0.01:+input0.prior;
  if(!(prior0>0&&prior0<0.5))throw new Error('swapseq: prior must be above 0 and below 0.5 (0.01 means 1%)');
  {const seen=new Set();for(const s of input0.samples){const id=String(s.id);if(seen.has(id))throw new Error(`swapseq: sample ID "${id}" appears more than once`);seen.add(id)}}
  // samples without a subject label are left out of the analysis and reported as not checked
  const missing=v=>v==null||['','na','nan','null','none','n/a','<na>'].includes(String(v).trim().toLowerCase());
  const keepIx=[];input0.subject.forEach((v,i)=>{if(!missing(v))keepIx.push(i)});
  const sub=a=>Array.isArray(a)?keepIx.map(i=>a[i]):a;
  const input={...input0,prior:prior0,samples:sub(input0.samples),subject:sub(input0.subject),house:sub(input0.house),well:sub(input0.well),plate:sub(input0.plate),host:sub(input0.host)};
  if(input.samples.length<4)throw new Error('swapseq: fewer than 4 samples have a subject label');
  // feature indices sorted and merged per sample (the distance code walks them in order)
  const tidy=s=>{const n=s.idx.length;let sorted=true;for(let k=1;k<n;k++)if(!(s.idx[k]>s.idx[k-1])){sorted=false;break}
    let I=s.idx,V=s.val;if(!sorted){const m=new Map();for(let k=0;k<n;k++)m.set(+s.idx[k],(m.get(+s.idx[k])||0)+(+s.val[k]));I=[...m.keys()].sort((x,y)=>x-y);V=I.map(i=>m.get(i))}
    const I2=[],V2=[];for(let k=0;k<I.length;k++){const v=+V[k];if(!Number.isFinite(v)||v<0)throw new Error(`swapseq: sample "${s.id}" has a negative or non-numeric value`);if(v>0){I2.push(+I[k]);V2.push(v)}}
    return {id:String(s.id),idx:Int32Array.from(I2),val:Float64Array.from(V2)}};
  const data={feats:input.features,samples:input.samples.map(tidy)};
  const lab=input.subject.map(v=>String(v).trim()),type=input.type==='markers'?'markers':'counts';
  const M=type==='markers'?markerDistances(data):distances(data);
  let res=annotate(data,lab,analyse(data,lab,{prior:prior0},M),{house:input.house||null,well:input.well||null,plate:input.plate||null});
  const mp=measurePower(data,lab,res,M); // v1.1: power measured on errors hidden in copies of the labels
  res.power=mp||{single:null,swap:null,singleWithCheck:null,swapWithCheck:null,measured:false};
  promotePairs(res,mp);secondOpinion(M,lab,res);
  wideScreen(data,lab,res,M,{house:input.house||null,well:input.well||null,plate:input.plate||null});
  const dup=findDuplicates(M,lab,{},data.samples.map(s=>s.id)).pairs; // two people's labels, one near-identical community
  let host=null;
  if(input.host){const G=input.host.map(h=>h?new Map(Object.entries(h).map(([k,v])=>[k,v])):null);
    // sites may be any string id; map them to integers once
    const ids=new Map();const G2=G.map(m=>{if(!m)return null;const o=new Map();for(const [k,v] of m){if(!ids.has(k))ids.set(k,ids.size);o.set(ids.get(k),v)}return o});
    host=hostCheck(G2,lab);res=addHost(res,host)}
  const num=v=>(v==null||Number.isNaN(v))?null:+v;
  const rows=res.rows.map((r,i)=>({sample:r.id,subject:lab[i],level:r.level,kind:r.kind||null,probability:num(r.prob),
    looks_like:r.looksLike==null?null:String(r.looksLike),partner:r.partner==null?null:String(r.partner),partner_subject:r.partnerSubject==null?null:String(r.partnerSubject),
    scorable:!!r.scorable,rare_total:r.rare?num(r.rare.total):null,rare_shared_own:r.rare?num(r.rare.withOwn):null,rare_shared_lookalike:r.rare?num(r.rare.withOther):null,household_lookalike:!!r.cohab,
    well:r.well==null?null:String(r.well),in_shift:!!r.shift,
    host_verdict:r.host?r.host.verdict:null,host_looks_like:r.host&&r.host.looksLike!=null?String(r.host.looksLike):null,host_llr:r.host?num(r.host.llr):null,host_sites:r.host?num(r.host.sites):null,
    flagged_by_host:!!r.byHost,flagged_as_pair:!!r.byPair,flagged_by_plate:!!r.byPlate,second_opinion:r.second==null?null:!!r.second,wide_screen:!!r.wide}));
  // put back the samples that had no subject label, in the caller's order
  let outRows=rows;
  if(keepIx.length!==input0.samples.length){const blank=id=>{const o={};for(const k in rows[0])o[k]=null;return Object.assign(o,{sample:String(id),level:'ok',scorable:false,household_lookalike:false,in_shift:false,flagged_by_host:false,flagged_as_pair:false,flagged_by_plate:false,second_opinion:null,wide_screen:false})};
    outRows=input0.samples.map(s=>blank(s.id));keepIx.forEach((i,k)=>outRows[i]=rows[k])}
  return JSON.stringify({rows:outRows,no_subject:input0.samples.length-keepIx.length,calibrated:!res.uncalibrated,shifts:res.shifts||[],duplicates:dup.map(p=>({sample_a:String(p.idA),subject_a:String(p.subjectA),sample_b:String(p.idB),subject_b:String(p.subjectB),distance:num(p.distance[0])})),power:res.power,n_subjects:res.nSubjects,separation:num(res.separation),promoted_pairs:res.promoted||0,swap_pairs:res.pairs.length,many_errors:res.nScorable>0&&res.pairs.length>=3&&2*res.pairs.length/res.nScorable>=1.5*prior0,separation2:num(res.separation2),reidentified:num(res.reid),n_scorable:res.nScorable,threshold:num(res.threshold),prior:prior0,
    type,markers:type==='markers'?M.markers:null,host:host?{ok:host.ok,reason:host.reason||null,samples_with_host_dna:host.used||0}:null,version:'1.1.1'});
}
