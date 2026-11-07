// Stress test of the host-DNA check under conditions real studies have (simulated genotypes; sites independent).
// Each scenario: people, samples, reads at common SNPs; hide label errors; count caught / named right / false flags.
import {hostCheck} from './core.mjs';
const SEED0=+(process.env.SEED||1);
function rng(seed){let x=seed>>>0||1;return()=>{x^=x<<13;x>>>=0;x^=x>>17;x^=x<<5;x>>>=0;return x/4294967296}}
function gauss(R){return Math.sqrt(-2*Math.log(R()+1e-12))*Math.cos(2*Math.PI*R())}

// genotypes for "people": each has {g: Uint8Array, pop}
function makeAF(R,K,lo=.05,hi=.5){return Float64Array.from({length:K},()=>lo+(hi-lo)*R())}
const draw=(R,af)=>Uint8Array.from(af,p=>(R()<p)+(R()<p));
function child(R,a,b){const K=a.length,c=new Uint8Array(K);for(let k=0;k<K;k++){const fa=a[k]===2?1:a[k]===0?0:(R()<.5?1:0),fb=b[k]===2?1:b[k]===0?0:(R()<.5?1:0);c[k]=fa+fb}return c}
// reads of a sample: nSites single reads at random panel sites; contamination c from a second genotype; read error e
function reads(R,K,g,nSites,e,cont){const m=new Map();for(let t=0;t<nSites;t++){const k=Math.floor(R()*K);const src=cont&&R()<cont.frac?cont.g:g;const v=src[k];let alt=v===2?1:v===0?0:(R()<.5?1:0);if(R()<e)alt=1-alt;m.set(k,[alt,1])}return m}

const ONLY=(process.env.ONLY||'').split(',').filter(Boolean);
function run(name,sc){if(ONLY.length&&!ONLY.includes(name.split(' ')[0]))return;
  const res={name,caught:0,named:0,swapped:0,ff:0,fN:0,unknown:0,n:0};
  for(let rep=0;rep<(sc.reps||5);rep++){const R=rng(SEED0*1000+rep*7919+name.length*31);
    const K=sc.K||200000,af=makeAF(R,K);
    // people
    const people=sc.people(R,K,af);              // [{id,g,fam?}]
    const G=[],lab=[],truth=[],owner=[];
    people.forEach((p,pi)=>{const per=typeof sc.per==='function'?sc.per(R,pi):(sc.per||6);
      for(let s=0;s<per;s++){const ns=sc.sites?sc.sites(R):3000;let cont=null;
        if(sc.cont&&R()<sc.cont.share){const other=people[(pi+1+Math.floor(R()*(people.length-1)))%people.length];cont={g:other.g,frac:sc.cont.frac}}
        G.push(reads(R,K,p.g,ns,sc.err??.01,cont));lab.push(p.id);truth.push(p.id);owner.push(pi)}});
    // hidden errors: swaps between random pairs, or within families when the scenario says so
    const L=lab.slice(),sw=new Set();const n=L.length;let guard=0;
    const pairs=sc.swapPairs?sc.swapPairs(R,people,owner):null;
    if(pairs)for(const [a,b] of pairs){[L[a],L[b]]=[L[b],L[a]];sw.add(a);sw.add(b)}
    else while(sw.size<(sc.nSwap||4)*2&&guard++<1e5){const a=Math.floor(R()*n),b=Math.floor(R()*n);if(L[a]===L[b]||sw.has(a)||sw.has(b))continue;[L[a],L[b]]=[L[b],L[a]];sw.add(a);sw.add(b)}
    const h=hostCheck(G,L,sc.opt||{});if(!h.ok){res.note=h.reason;continue}
    for(const r of h.rows){res.n++;if(r.verdict==='unknown')res.unknown++;const flag=r.verdict==='other person';if(r.verdict==='identical DNA')res.ident=(res.ident||0)+1;
      if(sw.has(r.i)){res.swapped++;if(flag){res.caught++;if(r.looksLike===truth[r.i])res.named++}}else{res.fN++;if(flag)res.ff++}}}
  const p=v=>res.swapped?Math.round(100*v/res.swapped)+'%':'n/a';
  console.log(`${name.padEnd(46)} caught ${p(res.caught).padStart(4)}  right person ${p(res.named).padStart(4)}  false flags/1000 ${(1000*res.ff/Math.max(1,res.fN)).toFixed(1).padStart(5)}  not judged ${Math.round(100*res.unknown/Math.max(1,res.n))}%  identical-DNA notes ${res.ident||0}${res.note?'  ('+res.note+')':''}`);
  return res}

const unrelated=(N)=>(R,K,af)=>Array.from({length:N},(_,i)=>({id:'P'+i,g:draw(R,af)}));
// families: mother + infant (+ father not sampled)
const motherInfant=(F)=>(R,K,af)=>{const out=[];for(let f=0;f<F;f++){const m=draw(R,af),d=draw(R,af),c=child(R,m,d);out.push({id:`M${f}`,g:m,fam:f},{id:`I${f}`,g:c,fam:f})}return out};
const siblings=(F)=>(R,K,af)=>{const out=[];for(let f=0;f<F;f++){const m=draw(R,af),d=draw(R,af);out.push({id:`S${f}a`,g:child(R,m,d),fam:f},{id:`S${f}b`,g:child(R,m,d),fam:f})}return out};
const twinsMZ=(F)=>(R,K,af)=>{const out=[];for(let f=0;f<F;f++){const g=draw(R,af);out.push({id:`T${f}a`,g,fam:f},{id:`T${f}b`,g:g.slice(),fam:f})}return out};
const twoPops=(N)=>(R,K,af)=>{const af2=Float64Array.from(af,p=>Math.min(.95,Math.max(.01,p+0.12*gauss(R))));return Array.from({length:N},(_,i)=>({id:'P'+i,g:draw(R,i<N/2?af:af2)}))};
// swap one sample of each of the first n families between the two relatives
const famSwap=n=>(R,people,owner)=>{const pairs=[];for(let f=0;f<n;f++){const a=owner.findIndex(o=>people[o].fam===f&&o%2===0),b=owner.findIndex(o=>people[o].fam===f&&o%2===1);if(a>=0&&b>=0)pairs.push([a,b])}return pairs};

console.log('Host-DNA check, stress test (simulated; 200,000-SNP panel; 5 repeats each)\n');
run('1 Baseline: 40 people x 6 samples, 3,000 sites',{people:unrelated(40)});
run('2 Uneven host DNA: sites lognormal, median 3,000',{people:unrelated(40),sites:R=>Math.round(3000*Math.exp(1.0*gauss(R)))});
run('3 Low host DNA: median 1,500 sites',{people:unrelated(40),sites:R=>Math.round(1500*Math.exp(.8*gauss(R)))});
run('4 Only 2 samples per person',{people:unrelated(60),per:2});
run('5 Only 3 samples per person',{people:unrelated(60),per:3});
run('6 Mixed: 2-8 samples per person',{people:unrelated(50),per:(R)=>2+Math.floor(R()*7)});
run('7 Mother-infant pairs, random swaps',{people:motherInfant(20)});
run('8 Mother-infant pairs, swaps WITHIN families',{people:motherInfant(20),swapPairs:famSwap(4)});
run('9 Siblings, swaps within families',{people:siblings(20),swapPairs:famSwap(4)});
run('10 Identical twins, swaps within twins',{people:twinsMZ(20),swapPairs:famSwap(4)});
run('11 Identical twins, random swaps (false flags?)',{people:twinsMZ(20)});
run('12 Two populations (allele frequencies differ)',{people:twoPops(40)});
run('13 Read error 3% (model assumes 1%)',{people:unrelated(40),err:.03});
run('14 Contamination: 10% of samples carry 20% foreign DNA',{people:unrelated(40),cont:{share:.1,frac:.2}});
run('15 Contamination: 10% of samples carry 45% foreign DNA',{people:unrelated(40),cont:{share:.1,frac:.45}});
run('16 Large cohort: 200 people x 4 samples',{people:unrelated(200),per:4,nSwap:10,reps:2});
run('17 Bigger panel (1M SNPs), same 3,000 sites',{people:unrelated(40),K:1000000});
