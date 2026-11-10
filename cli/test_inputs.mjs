// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mankirat Singh
// Robustness battery for the SwapSeq CLI (node cli/test_inputs.mjs): odd and broken input files. Every case: exit code, whether stderr is a clean one-line message or a stack trace,
// and for cases that should work, whether the planted swap (P1_1 <-> P2_1) is still the answer.
import {spawnSync} from 'child_process';import fs from 'fs';import os from 'os';import path from 'path';import {fileURLToPath} from 'url';
const CLI=fileURLToPath(new URL('./swapseq.mjs',import.meta.url)),D=fs.mkdtempSync(path.join(os.tmpdir(),'fz-'));
const R=(s=>()=>{s^=s<<13;s>>>=0;s^=s>>17;s^=s<<5;s>>>=0;return s/4294967296})(4242);for(let i=0;i<200;i++)R();const ri=(a,b)=>a+Math.floor(R()*(b-a+1));
function study(P=20,per=5,own=30,shared=50,idf=(p,k)=>`P${p}_${k}`,subf=p=>'P'+p){const F=shared+P*own,cols=[],ids=[],subj=[];
  for(let p=1;p<=P;p++)for(let k=1;k<=per;k++){ids.push(idf(p,k));subj.push(subf(p));const c=new Array(F).fill(0);
    for(let f=0;f<shared;f++)if(R()<.8)c[f]=ri(1,50);for(let f=0;f<own;f++)if(R()<.7)c[shared+(p-1)*own+f]=ri(1,30);cols.push(c)}
  return {F,cols,ids,subj,per}}
const S=study(30,5);const A_=10*S.per,B_=16*S.per;const sw=S.subj.slice();[sw[A_],sw[B_]]=[sw[B_],sw[A_]];
const table=(s,{sep='\t',eol='\n',head='#OTU ID',fmt=v=>v,bom='',trail='',pre='',q=x=>x}={})=>bom+pre+[head,...s.ids.map(q)].join(sep)+trail+eol+Array.from({length:s.F},(_,i)=>[q('f'+i),...s.cols.map(c=>fmt(c[i]))].join(sep)+trail).join(eol)+eol;
const meta=(s,lab,{sep='\t',eol='\n',head=['sample-id','person'],bom='',q=x=>x,extra=''}={})=>bom+head.join(sep)+eol+s.ids.map((id,i)=>[q(id),q(lab[i])].join(sep)).join(eol)+eol+extra;
const W=(n,t)=>{fs.writeFileSync(path.join(D,n),t);return n};
const cases=[];let small=null;const C=(name,expect,files,args)=>cases.push({name,expect,files,args});
// --- should work and find the swap ---
C('baseline','swap',()=>[W('t.tsv',table(S)),W('m.tsv',meta(S,sw))],['--subject','person']);
C('CRLF line ends','swap',()=>[W('t.tsv',table(S,{eol:'\r\n'})),W('m.tsv',meta(S,sw,{eol:'\r\n'}))],['--subject','person']);
C('old-Mac CR line ends','swap',()=>[W('t.tsv',table(S,{eol:'\r'})),W('m.tsv',meta(S,sw,{eol:'\r'}))],['--subject','person']);
C('UTF-8 BOM on both files','swap',()=>[W('t.tsv',table(S,{bom:'﻿'})),W('m.tsv',meta(S,sw,{bom:'﻿'}))],['--subject','person']);
C('CSV with every field quoted','swap',()=>[W('t.csv',table(S,{sep:',',q:x=>`"${x}"`})),W('m.csv',meta(S,sw,{sep:',',q:x=>`"${x}"`}))],['--subject','person']);
C('semicolon-separated (European Excel)','swap',()=>[W('t.csv',table(S,{sep:';'})),W('m.csv',meta(S,sw,{sep:';'}))],['--subject','person']);
C('R write.table header (no corner cell)','swap',()=>{const L=table(S).split('\n');L[0]=S.ids.join('\t');return [W('t.tsv',L.join('\n')),W('m.tsv',meta(S,sw))]},['--subject','person']);
C('tab-separated with every field quoted','swap',()=>[W('t.tsv',table(S,{q:x=>`"${x}"`})),W('m.tsv',meta(S,sw,{q:x=>`"${x}"`}))],['--subject','person']);
C('small study: 8 people x 4 samples, one swap','swap',()=>{const s2=study(8,4),l=s2.subj.slice();[l[8],l[20]]=[l[20],l[8]];small=[s2.ids[8],s2.ids[20]];return [W('t.tsv',table(s2)),W('m.tsv',meta(s2,l))]},['--subject','person']);
C('small study with clean labels: no flags','ok0',()=>{const s2=study(8,4);return [W('t.tsv',table(s2)),W('m.tsv',meta(s2,s2.subj))]},['--subject','person']);
C('trailing tab on every line','swap',()=>[W('t.tsv',table(S,{trail:'\t'})),W('m.tsv',meta(S,sw))],['--subject','person']);
C('biom-convert comment line first','swap',()=>[W('t.tsv',table(S,{pre:'# Constructed from biom file\n'})),W('m.tsv',meta(S,sw))],['--subject','person']);
C('blank lines at the end','swap',()=>[W('t.tsv',table(S)+'\n\n\n'),W('m.tsv',meta(S,sw)+'\n\n')],['--subject','person']);
C('decimal counts 12.0','swap',()=>[W('t.tsv',table(S,{fmt:v=>v.toFixed(1)})),W('m.tsv',meta(S,sw))],['--subject','person']);
C('relative abundances (0-1)','swap',()=>{const tot=S.cols.map(c=>c.reduce((a,b)=>a+b,0));const s2={...S,cols:S.cols.map((c,j)=>c.map(v=>v/tot[j]))};return [W('t.tsv',table(s2,{fmt:v=>v?v.toPrecision(6):'0'})),W('m.tsv',meta(S,sw))]},['--subject','person']);
C('scientific notation 1.2e+01','swap',()=>[W('t.tsv',table(S,{fmt:v=>v.toExponential(2)})),W('m.tsv',meta(S,sw))],['--subject','person']);
C('empty cells instead of 0','swap',()=>[W('t.tsv',table(S,{fmt:v=>v||''})),W('m.tsv',meta(S,sw))],['--subject','person']);
C('taxonomy column at the end','swap',()=>{const t=table(S).split('\n').filter(Boolean).map((l,i)=>l+'\t'+(i?'k__Bacteria; p__Firmicutes':'taxonomy')).join('\n')+'\n';return [W('t.tsv',t),W('m.tsv',meta(S,sw))]},['--subject','person']);
C('metadata rows in another order','swap',()=>{const o=S.ids.map((_,i)=>i).reverse();return [W('t.tsv',table(S)),W('m.tsv','sample-id\tperson\n'+o.map(i=>S.ids[i]+'\t'+sw[i]).join('\n')+'\n')]},['--subject','person']);
C('metadata has extra samples not in table','swap',()=>[W('t.tsv',table(S)),W('m.tsv',meta(S,sw,{extra:'ZZ_1\tP99\nZZ_2\tP99\n'}))],['--subject','person']);
C('QIIME 2 #q2:types row in metadata','swap',()=>{const L=meta(S,sw).split('\n');L.splice(1,0,'#q2:types\tcategorical');return [W('t.tsv',table(S)),W('m.tsv',L.join('\n'))]},['--subject','person']);
C('spaces around IDs and labels','swap',()=>[W('t.tsv',table(S,{q:x=>' '+x+' '})),W('m.tsv',meta(S,sw,{q:x=>x+' '}))],['--subject','person']);
C('numeric subject labels (1,2,3)','swap',()=>[W('t.tsv',table(S)),W('m.tsv',meta(S,sw.map(x=>x.slice(1))))],['--subject','person']);
C('numeric sample IDs','swap',()=>{const s2={...S,ids:S.ids.map((_,i)=>String(1000+i))};return [W('t.tsv',table(s2)),W('m.tsv',meta(s2,sw))]},['--subject','person']);
C('IDs with dots, dashes, unicode','swap',()=>{const s2={...S,ids:S.ids.map((x,i)=>`Patïent-${x}.v1`)};return [W('t.tsv',table(s2)),W('m.tsv',meta(s2,sw.map(x=>'José '+x)))]},['--subject','person']);
C('subject column name with space, case differs','swap',()=>[W('t.tsv',table(S)),W('m.tsv',meta(S,sw,{head:['#SampleID','Subject ID']}))],['--subject','subject id']);
C('samples in rows','swap',()=>{const t='sample\t'+Array.from({length:S.F},(_,i)=>'f'+i).join('\t')+'\n'+S.ids.map((id,j)=>id+'\t'+S.cols[j].join('\t')).join('\n')+'\n';return [W('t.tsv',t),W('m.tsv',meta(S,sw))]},['--subject','person','--samples-in-rows']);
C('samples in rows WITHOUT the flag (auto-detect or clear error)','swap|err',()=>{const t='sample\t'+Array.from({length:S.F},(_,i)=>'f'+i).join('\t')+'\n'+S.ids.map((id,j)=>id+'\t'+S.cols[j].join('\t')).join('\n')+'\n';return [W('t.tsv',t),W('m.tsv',meta(S,sw))]},['--subject','person']);
C('all-zero feature rows and one all-zero sample','swap',()=>{const s2={...S,cols:S.cols.map((c,j)=>j===120?c.map(()=>0):c)};return [W('t.tsv',table(s2)),W('m.tsv',meta(S,sw))]},['--subject','person']);
C('some people have only 1 sample','swap',()=>{const lab=sw.slice();for(let p=20;p<30;p++)for(let k=1;k<S.per;k++)lab[p*S.per+k]='solo'+p+'_'+k;return [W('t.tsv',table(S)),W('m.tsv',meta(S,lab))]},['--subject','person']);
C('missing subject for some samples (NA / empty)','swap',()=>{const lab=sw.slice();lab[100]='NA';lab[101]='';lab[102]='nan';return [W('t.tsv',table(S)),W('m.tsv',meta(S,lab))]},['--subject','person']);
C('BIOM 1.0 JSON (sparse)','swap',()=>{const data=[];S.cols.forEach((c,j)=>c.forEach((v,i)=>{if(v)data.push([i,j,v])}));return [W('t.biom',JSON.stringify({format:'Biological Observation Matrix 1.0.0',matrix_type:'sparse',shape:[S.F,S.ids.length],rows:Array.from({length:S.F},(_,i)=>({id:'f'+i})),columns:S.ids.map(id=>({id})),data})),W('m.tsv',meta(S,sw))]},['--subject','person']);
C('BIOM 1.0 JSON (dense)','swap',()=>[W('t.biom',JSON.stringify({format:'Biological Observation Matrix 1.0.0',matrix_type:'dense',shape:[S.F,S.ids.length],rows:Array.from({length:S.F},(_,i)=>({id:'f'+i})),columns:S.ids.map(id=>({id})),data:Array.from({length:S.F},(_,i)=>S.cols.map(c=>c[i]))})),W('m.tsv',meta(S,sw))],['--subject','person']);
// --- should stop with a clear message (exit 1, one or two lines, no stack trace) ---
C('duplicate sample ID in table','err',()=>{const s2={...S,ids:S.ids.map((x,i)=>i===7?S.ids[6]:x)};return [W('t.tsv',table(s2)),W('m.tsv',meta(S,sw))]},['--subject','person']);
C('duplicate sample ID in metadata with DIFFERENT subjects','err',()=>[W('t.tsv',table(S)),W('m.tsv',meta(S,sw,{extra:'P3_2\tP17\n'}))],['--subject','person']);
C('text where a number should be','err',()=>[W('t.tsv',table(S,{fmt:v=>v===7?'seven':v})),W('m.tsv',meta(S,sw))],['--subject','person']);
C('negative counts','err',()=>[W('t.tsv',table(S,{fmt:v=>v===7?-7:v})),W('m.tsv',meta(S,sw))],['--subject','person']);
C('NaN / Inf values','err',()=>[W('t.tsv',table(S,{fmt:v=>v===7?'NaN':v===8?'Inf':v})),W('m.tsv',meta(S,sw))],['--subject','person']);
C('European decimal comma 0,5','err|swap',()=>[W('t.tsv',table(S,{fmt:v=>String(v+.5).replace('.',',')})),W('m.tsv',meta(S,sw))],['--subject','person']);
C('no sample ID overlaps','err',()=>{const s2={...S,ids:S.ids.map(x=>'X'+x)};return [W('t.tsv',table(s2)),W('m.tsv',meta(S,sw))]},['--subject','person']);
C('only the first third of samples have metadata (swap not among them)','ok',()=>[W('t.tsv',table(S)),W('m.tsv','sample-id\tperson\n'+S.ids.slice(0,50).map((id,i)=>id+'\t'+sw[i]).join('\n')+'\n')],['--subject','person']);
C('subject column missing','err',()=>[W('t.tsv',table(S)),W('m.tsv',meta(S,sw))],['--subject','patient']);
C('every sample a different person','err',()=>[W('t.tsv',table(S)),W('m.tsv',meta(S,S.ids))],['--subject','person']);
C('everyone the same person','err',()=>[W('t.tsv',table(S)),W('m.tsv',meta(S,S.ids.map(()=>'A')))],['--subject','person']);
C('only 2 people','err|ok',()=>{const s2=study(2,6);return [W('t.tsv',table(s2)),W('m.tsv',meta(s2,s2.subj))]},['--subject','person']);
C('empty table file','err',()=>[W('t.tsv',''),W('m.tsv',meta(S,sw))],['--subject','person']);
C('header only, no features','err',()=>[W('t.tsv','#OTU ID\t'+S.ids.join('\t')+'\n'),W('m.tsv',meta(S,sw))],['--subject','person']);
C('table file missing','err',()=>['nope.tsv',W('m.tsv',meta(S,sw))],['--subject','person']);
C('table and metadata swapped on the command line','err',()=>{W('t.tsv',table(S));W('m.tsv',meta(S,sw));return ['m.tsv','t.tsv']},['--subject','person']);
C('ragged rows (a line with fewer columns)','err|swap',()=>{const L=table(S).split('\n');L[5]=L[5].split('\t').slice(0,40).join('\t');return [W('t.tsv',L.join('\n')),W('m.tsv',meta(S,sw))]},['--subject','person']);
C('binary garbage as table','err',()=>{fs.writeFileSync(path.join(D,'t.tsv'),Buffer.from(Array.from({length:5000},()=>ri(0,255))));return ['t.tsv',W('m.tsv',meta(S,sw))]},['--subject','person']);
C('an Excel .xlsx given as table','err',()=>{fs.writeFileSync(path.join(D,'t.xlsx'),Buffer.concat([Buffer.from('PK\u0003\u0004'),Buffer.from(Array.from({length:3000},()=>ri(0,255)))]));return ['t.xlsx',W('m.tsv',meta(S,sw))]},['--subject','person']);
C('BIOM 2 (HDF5) given as table','err',()=>{fs.writeFileSync(path.join(D,'t.biom'),Buffer.concat([Buffer.from([0x89,0x48,0x44,0x46,0x0d,0x0a,0x1a,0x0a]),Buffer.from(Array.from({length:3000},()=>ri(0,255)))]));return ['t.biom',W('m.tsv',meta(S,sw))]},['--subject','person']);
C('--prior 0','err',()=>[W('t.tsv',table(S)),W('m.tsv',meta(S,sw))],['--subject','person','--prior','0']);
C('--prior 5 (meant 5%)','err',()=>[W('t.tsv',table(S)),W('m.tsv',meta(S,sw))],['--subject','person','--prior','5']);
C('--prior abc','err',()=>[W('t.tsv',table(S)),W('m.tsv',meta(S,sw))],['--subject','person','--prior','abc']);
C('--min-depth removes everything','err',()=>[W('t.tsv',table(S)),W('m.tsv',meta(S,sw))],['--subject','person','--min-depth','1000000']);
C('--type wrong','err',()=>[W('t.tsv',table(S)),W('m.tsv',meta(S,sw))],['--subject','person','--type','otu']);
C('--out would overwrite the input table','err',()=>[W('t.tsv',table(S)),W('m.tsv',meta(S,sw))],['--subject','person','--out','t']);
C('--out into a folder that does not exist','err',()=>[W('t.tsv',table(S)),W('m.tsv',meta(S,sw))],['--subject','person','--out','nodir/x/res']);
// run
let bad=0;const only=process.argv[2];
for(const c of cases){if(only&&!c.name.includes(only))continue;for(const f of fs.readdirSync(D))fs.rmSync(path.join(D,f),{recursive:true,force:true});
  const files=c.files();const out=c.args.includes('--out')?[]:['--out','r'];
  const r=spawnSync('node',[CLI,...files,...c.args,...out,'--no-power'],{cwd:D,encoding:'utf8',timeout:120000});
  const trace=/\n\s+at .*\(|node:internal|TypeError|RangeError|ReferenceError|SyntaxError|Maximum call/.test(r.stderr||'');
  let got;if(r.error)got='TIMEOUT';else if(r.status===0||r.status===2){let fl=null;try{fl=fs.readFileSync(path.join(D,'r.tsv'),'utf8').split('\n').slice(2).filter(Boolean).map(l=>l.split('\t')).filter(x=>x[3]==='flag').map(x=>x[0].trim())}catch(e){}
      const nm=i=>c.name.includes('numeric sample IDs')?String(1000+i):c.name.includes('unicode')?`Pat\u00efent-${S.ids[i]}.v1`:S.ids[i];const want=c.name.startsWith('small study: ')?small:[nm(A_),nm(B_)];
      got=fl&&fl.length===2&&want.every(w=>fl.includes(w))?'swap':fl&&fl.length===0&&c.expect==='ok0'?'ok0':'ok('+(fl?fl.length+' flags: '+fl.slice(0,4).join(','):'no output')+')'}
    else got=trace?'CRASH':'err';
  const pass=c.expect.split('|').some(e=>e===got||(e==='ok'&&got.startsWith('ok')));if(!pass)bad++;
  console.log(`${pass?'  ok ':'FAIL '} ${c.name.padEnd(58)} want ${c.expect.padEnd(9)} got ${got}${got==='err'||got==='CRASH'||!pass?'  | '+(r.stderr||r.stdout||'').split('\n').filter(l=>l&&!l.startsWith('Run swapseq')).slice(0,2).join(' / ').slice(0,170):''}`)}
console.log(`\n${cases.length-bad} of ${cases.length} as expected`);process.exit(bad?1:0);
