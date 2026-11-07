// SPDX-License-Identifier: GPL-3.0-or-later
// Tests for the SwapSeq command-line tool on small simulated studies. Run: node cli/test.mjs
import {execFileSync} from 'child_process';import fs from 'fs';import os from 'os';import path from 'path';import assert from 'assert';
const CLI=new URL('./swapseq.mjs',import.meta.url).pathname,D=fs.mkdtempSync(path.join(os.tmpdir(),'swapseq-'));
const R=(s=>()=>{s^=s<<13;s>>>=0;s^=s>>17;s^=s<<5;s>>>=0;return s/4294967296})(12345);
const ri=(a,b)=>a+Math.floor(R()*(b-a+1));
// 20 people x 5 samples: 50 shared taxa + 40 personal taxa each
const P=20,per=5,own=40,shared=50,F=shared+P*own,ids=[],subj=[],cols=[];
for(let p=1;p<=P;p++)for(let k=1;k<=per;k++){ids.push(`P${p}_${k}`);subj.push('P'+p);const c=new Array(F).fill(0);
  for(let f=0;f<shared;f++)if(R()<.8)c[f]=ri(1,50);for(let f=0;f<own;f++)if(R()<.7)c[shared+(p-1)*own+f]=ri(1,30);cols.push(c)}
const swapped=subj.slice();[swapped[0],swapped[per]]=[swapped[per],swapped[0]];   // P1_1 <-> P2_1
const tsv=f=>fs.writeFileSync(f,'#OTU ID\t'+ids.join('\t')+'\n'+Array.from({length:F},(_,i)=>['f'+i,...cols.map(c=>c[i])].join('\t')).join('\n')+'\n');
const meta=(f,lab,sep='\t')=>fs.writeFileSync(f,['sample-id','person'].join(sep)+'\n'+ids.map((id,i)=>[id,lab[i]].join(sep)).join('\n')+'\n');
const run=(args,expectCode=0)=>{try{const o=execFileSync('node',[CLI,...args,'--quiet'],{cwd:D,encoding:'utf8'});assert.equal(0,expectCode,'expected exit '+expectCode);return o}catch(e){if(e.status!==expectCode)throw new Error(`exit ${e.status}: ${e.stderr||e.message}`);return e.stdout}};
const flagged=prefix=>fs.readFileSync(path.join(D,prefix+'.tsv'),'utf8').split('\n').slice(2).filter(Boolean).map(l=>l.split('\t')).filter(c=>c[3]==='flag').map(c=>c[0]).sort();
let ok=0;const test=(name,fn)=>{fn();ok++;console.log('ok',name)};

test('TSV table: finds the planted swap and nothing else',()=>{tsv(path.join(D,'t.tsv'));meta(path.join(D,'m.tsv'),swapped);
  run(['t.tsv','m.tsv','--subject','person','--out','a']);assert.deepEqual(flagged('a'),['P1_1','P2_1']);
  const j=JSON.parse(fs.readFileSync(path.join(D,'a.json'),'utf8'));assert.equal(j.strongFlags,2);assert.ok(j.power.measured);assert.match(j.methodsSentence,/SwapSeq 1\.1\.1/)});
test('clean labels: no strong flags',()=>{meta(path.join(D,'clean.tsv'),subj);run(['t.tsv','clean.tsv','--subject','person','--out','b','--no-power']);assert.deepEqual(flagged('b'),[])});
test('BIOM 1.0 JSON and CSV metadata give the same answer',()=>{const data=[];cols.forEach((c,j)=>c.forEach((v,i)=>{if(v)data.push([i,j,v])}));
  fs.writeFileSync(path.join(D,'t.biom'),JSON.stringify({id:'x',format:'Biological Observation Matrix 1.0.0',type:'OTU table',matrix_type:'sparse',shape:[F,ids.length],rows:Array.from({length:F},(_,i)=>({id:'f'+i,metadata:null})),columns:ids.map(id=>({id,metadata:null})),data}));
  meta(path.join(D,'m.csv'),swapped,',');run(['t.biom','m.csv','--subject','person','--out','c','--no-power']);assert.deepEqual(flagged('c'),['P1_1','P2_1'])});
test('--fail-on-flags exits 2',()=>{run(['t.tsv','m.tsv','--subject','person','--out','d','--no-power','--fail-on-flags'],2)});
test('markers mode',()=>{
  // strain-marker tables have many markers: 300 shared + 150 per person, each present or absent
  const MF=300+P*150,mk=ids.map((_,j)=>{const p=+subj[j].slice(1),c=new Uint8Array(MF);for(let f=0;f<300;f++)c[f]=R()<.5;for(let f=0;f<150;f++)c[300+(p-1)*150+f]=R()<.7;return c});
  fs.writeFileSync(path.join(D,'mk.tsv'),'marker\t'+ids.join('\t')+'\n'+Array.from({length:MF},(_,i)=>['m'+i,...mk.map(c=>c[i])].join('\t')).join('\n')+'\n');
  run(['mk.tsv','m.tsv','--subject','person','--type','markers','--out','e','--no-power']);assert.deepEqual(flagged('e'),['P1_1','P2_1'])});
test('host DNA catches a swap the microbes cannot see',()=>{
  // microbes carry no personal signal at all; a second swap (P3_1 <-> P4_1) is only visible in host DNA
  const flat=ids.map((_,j)=>cols[j].map((v,i)=>i<shared?v:0));
  fs.writeFileSync(path.join(D,'flat.tsv'),'#OTU ID\t'+ids.join('\t')+'\n'+Array.from({length:shared},(_,i)=>['f'+i,...flat.map(c=>c[i])].join('\t')).join('\n')+'\n');
  const lab=subj.slice();[lab[2*per],lab[3*per]]=[lab[3*per],lab[2*per]];meta(path.join(D,'mh.tsv'),lab);
  const K=200000,af=Float64Array.from({length:K},()=>.05+.45*R()),geno={};for(let p=1;p<=P;p++)geno['P'+p]=Uint8Array.from(af,q=>(R()<q)+(R()<q));
  const lines=['sample\tsite\talt\ttotal'];ids.forEach((id,j)=>{const seen=new Set();for(let t=0;t<4000;t++){const k=Math.floor(R()*K);if(seen.has(k))continue;seen.add(k);const g=geno[subj[j]][k];lines.push(`${id}\trs${k}\t${g===2?1:g===0?0:(R()<.5?1:0)}\t1`)}});
  fs.writeFileSync(path.join(D,'host.tsv'),lines.join('\n')+'\n');
  run(['flat.tsv','mh.tsv','--subject','person','--host','host.tsv','--out','h','--no-power']);assert.deepEqual(flagged('h'),['P3_1','P4_1']);
  const j=JSON.parse(fs.readFileSync(path.join(D,'h.json'),'utf8'));assert.ok(j.host.ok);assert.equal(j.host.flaggedByHost,2)});
test('helpful errors',()=>{run(['t.tsv','m.tsv','--subject','nope'],1);run(['missing.tsv','m.tsv','--subject','person'],1);run(['t.tsv'],1)});
console.log(`\n${ok} tests passed`);fs.rmSync(D,{recursive:true,force:true});
