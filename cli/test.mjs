// SPDX-License-Identifier: GPL-3.0-or-later
// Tests for the SwapSeq command-line tool on small simulated studies. Run: node cli/test.mjs
import {execFileSync} from 'child_process';import fs from 'fs';import os from 'os';import path from 'path';import assert from 'assert';import {fileURLToPath} from 'url';
const CLI=fileURLToPath(new URL('./swapseq.mjs',import.meta.url)),D=fs.mkdtempSync(path.join(os.tmpdir(),'swapseq-'));
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
  const j=JSON.parse(fs.readFileSync(path.join(D,'a.json'),'utf8'));assert.equal(j.strongFlags,2);assert.ok(j.power.measured);assert.match(j.methodsSentence,/SwapSeq 1\.1\.1/);
  const h=fs.readFileSync(path.join(D,'a.html'),'utf8');assert.match(h,/Sample identity report/);assert.match(h,/P1_1/);assert.ok(!/<script/i.test(h),'report has no scripts')});
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
test('duplicate check: one tube under two names',()=>{
  // P5_3 is replaced by a re-sequenced copy of P6_2 (every count jittered a little), still labelled P5
  const dup=cols.map(c=>c.slice()),src=ids.indexOf('P6_2'),dst=ids.indexOf('P5_3');dup[dst]=cols[src].map(v=>v?Math.max(1,Math.round(v*(0.85+0.3*R()))):0);
  fs.writeFileSync(path.join(D,'dup.tsv'),'#OTU ID\t'+ids.join('\t')+'\n'+Array.from({length:F},(_,i)=>['f'+i,...dup.map(c=>c[i])].join('\t')).join('\n')+'\n');
  run(['dup.tsv','clean.tsv','--subject','person','--out','u','--no-power']);
  const j=JSON.parse(fs.readFileSync(path.join(D,'u.json'),'utf8'));
  assert.deepEqual(j.duplicates.map(d=>[d.sampleA,d.sampleB].sort()),[['P5_3','P6_2']]);
  const b=JSON.parse(fs.readFileSync(path.join(D,'b.json'),'utf8'));assert.equal(b.duplicates.length,0)});
test('plate map, second check and wide screen leave the main result alone',()=>{
  const wells=ids.map((_,i)=>String.fromCharCode(65+i%8)+String(Math.floor(i/8)%12+1).padStart(2,'0')),plates=ids.map((_,i)=>'p'+(1+Math.floor(i/96)));
  // put the two swapped tubes side by side: P1_1 in A01 (already), P2_1 moved to B01
  const j=ids.indexOf('P2_1'),k=wells.indexOf('B01');[wells[j],wells[k]]=[wells[k],wells[j]];
  fs.writeFileSync(path.join(D,'mw.tsv'),'sample-id\tperson\twell\tplate\n'+ids.map((id,i)=>[id,swapped[i],wells[i],plates[i]].join('\t')).join('\n')+'\n');
  run(['t.tsv','mw.tsv','--subject','person','--well','well','--plate','plate','--out','w','--no-power']);
  assert.deepEqual(flagged('w'),['P1_1','P2_1']);
  const L=fs.readFileSync(path.join(D,'w.tsv'),'utf8').split('\n').filter(Boolean),h=L[1].split('\t'),rows=L.slice(2).map(l=>l.split('\t'));
  for(const c of ['flagged_by_plate','second_opinion','wide_screen'])assert.ok(h.includes(c),'column '+c);
  const col=n=>h.indexOf(n);
  for(const r of rows){if(r[col('wide_screen')]==='yes')assert.equal(r[col('level')],'ok','a wide-screen sample is not on the main lists');
    if(r[col('second_opinion')])assert.notEqual(r[col('level')],'ok','the second check only annotates listed samples')}
  for(const id of ['P1_1','P2_1'])assert.equal(rows.find(r=>r[0]===id)[col('second_opinion')],'agrees');
  const s=JSON.parse(fs.readFileSync(path.join(D,'w.json'),'utf8'));assert.ok(s.plate&&s.plate.used,'plate map used');assert.ok(s.wideScreen&&Array.isArray(s.wideScreen.extraSamples));
  // two samples in one well: the plate map is not trusted, and the result is the same as without it
  fs.writeFileSync(path.join(D,'mw2.tsv'),'sample-id\tperson\twell\n'+ids.map((id,i)=>[id,swapped[i],'A01'].join('\t')).join('\n')+'\n');
  run(['t.tsv','mw2.tsv','--subject','person','--well','well','--out','w2','--no-power','--no-wide']);assert.deepEqual(flagged('w2'),['P1_1','P2_1']);
  const s2=JSON.parse(fs.readFileSync(path.join(D,'w2.json'),'utf8'));assert.equal(s2.plate.used,false);assert.equal(s2.wideScreen,null)});
test('helpful errors',()=>{run(['t.tsv','m.tsv','--subject','nope'],1);run(['missing.tsv','m.tsv','--subject','person'],1);run(['t.tsv'],1)});
console.log(`\n${ok} tests passed`);fs.rmSync(D,{recursive:true,force:true});
