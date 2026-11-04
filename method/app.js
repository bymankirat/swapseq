// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mankirat Singh
// ---------- state ----------
const $=id=>document.getElementById(id);
const st={table:null,meta:null,tableName:'',metaName:''};
const esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmtPct=v=>Number.isFinite(v)?Math.round(v*100)+'%':'–';
function setStatus(msg,err){$('status').innerHTML=err?`<span class="err">${esc(msg)}</span>`:esc(msg)}

// ---------- parsing ----------
// A table is either {kind:'sparse', feats, samples:[{id,idx,val}]} or {kind:'grid', rowIds, colIds, get(r,c)} (orientation decided later)
async function readTable(file){
  const buf=new Uint8Array(await file.arrayBuffer());return parseTableBytes(buf,file.name)}
async function parseTableBytes(buf,name){
  if(buf[0]===0x50&&buf[1]===0x4b){ // .qza is a zip
    const zip=await JSZip.loadAsync(buf);const files=Object.keys(zip.files);
    const pick=files.find(f=>/\/data\/[^/]+\.biom$/.test(f))||files.find(f=>/\/data\/[^/]+\.(tsv|txt|csv)$/.test(f));
    if(!pick)throw new Error('This .qza has no feature table inside (looked for data/*.biom).');
    const typ=files.find(f=>/metadata\.yaml$/.test(f)&&f.split('/').length===2);
    if(typ){const y=await zip.file(typ).async('string');const m=y.match(/type:\s*(.+)/);if(m&&!/FeatureTable\[(Frequency|RelativeFrequency)\]/.test(m[1]))throw new Error(`This artifact is ${m[1].trim()}, not a FeatureTable[Frequency].`)}
    return parseTableBytes(await zip.file(pick).async('uint8array'),pick)}
  if(buf[0]===0x89&&buf[1]===0x48&&buf[2]===0x44&&buf[3]===0x46){ // HDF5 = BIOM 2.x
    const f=new hdf5.File(buf.buffer.slice(buf.byteOffset,buf.byteOffset+buf.byteLength),name);
    const g=p=>f.get(p).value;const feats=g('observation/ids'),ids=g('sample/ids'),ip=g('sample/matrix/indptr'),ix=g('sample/matrix/indices'),dt=g('sample/matrix/data');
    const samples=ids.map((id,j)=>{const a=ip[j],b=ip[j+1],o=[];for(let k=a;k<b;k++)if(dt[k]>0)o.push([ix[k],dt[k]]);o.sort((x,y)=>x[0]-y[0]);return {id:String(id),idx:o.map(x=>x[0]),val:o.map(x=>x[1])}});
    return {kind:'sparse',feats:feats.map(String),samples}}
  const txt=new TextDecoder().decode(buf);const t0=txt.trimStart();
  if(t0[0]==='{'){ // BIOM 1.0 JSON
    const j=JSON.parse(t0);const feats=j.rows.map(r=>r.id),ids=j.columns.map(c=>c.id),per=ids.map(()=>[]);
    if(j.matrix_type==='dense'){j.data.forEach((row,r)=>row.forEach((v,c)=>{if(v>0)per[c].push([r,v])}))}else for(const [r,c,v] of j.data)if(v>0)per[c].push([r,v]);
    return {kind:'sparse',feats,samples:ids.map((id,c)=>{const a=per[c].sort((x,y)=>x[0]-y[0]);return {id,idx:a.map(x=>x[0]),val:a.map(x=>x[1])}})}}
  // delimited text
  let lines=txt.split(/\r?\n/).filter(l=>l.trim()!=='');
  while(lines.length&&/^#\s*Constructed from biom/i.test(lines[0]))lines.shift();
  const delim=(lines[0].match(/\t/g)||[]).length>=(lines[0].match(/,/g)||[]).length?'\t':',';
  const split=l=>delim==='\t'?l.split('\t'):splitCSV(l);
  const head=split(lines[0]).map(s=>s.trim());const body=lines.slice(1).map(split);
  // drop non-numeric columns such as taxonomy
  const numericCol=head.map((_,c)=>c===0?false:body.slice(0,50).every(r=>r[c]===undefined||r[c]===''||isFinite(+r[c])));
  const cols=head.map((h,c)=>c).filter(c=>numericCol[c]);
  if(!cols.length)throw new Error('Could not find numeric count columns in this table.');
  const rowIds=body.map(r=>String(r[0]).trim()),colIds=cols.map(c=>head[c]);
  const M=body.map(r=>cols.map(c=>+r[c]||0));
  return {kind:'grid',rowIds,colIds,get:(r,c)=>M[r][c]}}
function splitCSV(l){const o=[];let cur='',q=false;for(let i=0;i<l.length;i++){const ch=l[i];if(q){if(ch==='"'){if(l[i+1]==='"'){cur+='"';i++}else q=false}else cur+=ch}else if(ch==='"')q=true;else if(ch===','){o.push(cur);cur=''}else cur+=ch}o.push(cur);return o}

async function readMeta(file){return parseMeta(await file.text())}
function parseMeta(txt){
  let lines=txt.split(/\r?\n/).filter(l=>l.trim()!=='');
  const delim=(lines[0].match(/\t/g)||[]).length>=(lines[0].match(/,/g)||[]).length?'\t':',';
  const split=l=>(delim==='\t'?l.split('\t'):splitCSV(l)).map(s=>s.trim());
  const head=split(lines[0]).map(h=>h.replace(/^#/,''));lines=lines.slice(1).filter(l=>!/^#q2:types/i.test(l));
  const rows=lines.filter(l=>!l.startsWith('#')).map(split);
  const map=new Map();rows.forEach(r=>map.set(r[0],r));
  return {head,rows,map}}

// orient a table against metadata and make the sparse sample list
function samplesFrom(table,meta){
  if(table.kind==='sparse')return {feats:table.feats.length,samples:table.samples};
  const colHit=table.colIds.filter(id=>meta.map.has(id)).length,rowHit=table.rowIds.filter(id=>meta.map.has(id)).length;
  if(colHit>=rowHit){ // features in rows (QIIME/BIOM default)
    return {feats:table.rowIds.length,samples:table.colIds.map((id,c)=>{const idx=[],val=[];for(let r=0;r<table.rowIds.length;r++){const v=table.get(r,c);if(v>0){idx.push(r);val.push(v)}}return {id,idx,val}})}}
  return {feats:table.colIds.length,samples:table.rowIds.map((id,r)=>{const idx=[],val=[];for(let c=0;c<table.colIds.length;c++){const v=table.get(r,c);if(v>0){idx.push(c);val.push(v)}}return {id,idx,val}})}
}

// ---------- column guessing ----------
const PERSON=/^(host_?subject_?id|subject(_?id)?|participant(_?id)?|patient(_?id)?|personal_?id|individual(_?id)?|host_?id|person(_?id)?|donor(_?id)?|studyid|study_?id|mouse(_?id)?|animal(_?id)?|child(_?id)?|volunteer(_?id)?|pid)$/i;
const GROUP=/^(body_?site|bodysite|sample_?type|site|habitat|body_?habitat|env_?material|tissue|sample_?site)$/i;
function fillColumns(){
  const m=st.meta,P=$('colP'),G=$('colG');const cols=m.head.slice(1);
  const nVals=c=>{const k=m.head.indexOf(c);return new Set(m.rows.map(r=>r[k])).size};
  P.innerHTML=cols.map(c=>`<option>${esc(c)}</option>`).join('');G.innerHTML='<option value="">No split</option>'+cols.map(c=>`<option>${esc(c)}</option>`).join('');
  const guessP=cols.find(c=>PERSON.test(c))||cols.find(c=>{const v=nVals(c);return v>1&&v<m.rows.length*.8&&/id|subj|person|part|pat/i.test(c)})||cols[0];
  P.value=guessP;const guessG=cols.find(c=>GROUP.test(c)&&nVals(c)>1&&nVals(c)<=12);G.value=guessG||'';
  const opt=cols.map(c=>`<option>${esc(c)}</option>`).join('');
  for(const [id,lbl] of [['colH','None'],['colW','None'],['colL','None']]){$(id).innerHTML=`<option value="">${lbl}</option>`+opt;$(id).disabled=false}
  const vals=c=>{const k=m.head.indexOf(c);return m.rows.map(r=>r[k]).filter(Boolean)};
  $('colH').value=cols.find(c=>/^(household|house(hold)?_?id|family(_?id)?|home(_?id)?|cohab\w*|room(mate)?s?(_?id)?)$/i.test(c))||'';
  $('colW').value=cols.find(c=>{const v=vals(c);return v.length&&v.filter(x=>/^[A-P]\s*0?\d{1,2}$/i.test(x.trim())).length>=.9*v.length&&/well|pos/i.test(c)})||cols.find(c=>{const v=vals(c);return v.length>8&&v.filter(x=>/^[A-H]0?\d{1,2}$/i.test(x.trim())).length>=.95*v.length})||'';
  $('colL').value=$('colW').value?(cols.find(c=>/plate/i.test(c)&&c!==$('colW').value)||''):'';
  P.disabled=G.disabled=false;document.querySelectorAll('.colfield').forEach(e=>e.hidden=false);
}
function stat(id,ok,txt){const li=$(id);if(!li)return;li.classList.toggle('ok',ok);li.querySelector('.st').textContent=txt}
function ready(){
  document.querySelectorAll('.colfield').forEach(e=>e.hidden=!st.meta);
  if(st.table&&st.meta){const s=samplesFrom(st.table,st.meta);const hit=s.samples.filter(x=>st.meta.map.has(x.id)).length;
    $('infoM').innerHTML=`${esc(st.metaName)}: <b>${hit}</b> of ${s.samples.length} table samples found in the metadata.`;
    stat('statM',hit>0,hit?`${hit.toLocaleString()} of ${s.samples.length.toLocaleString()} samples matched`:'no sample IDs match the table');
    if(!hit)setStatus('No sample IDs match between the table and the metadata. Check the first column of the metadata.',true);else setStatus('')}
}

// ---------- file inputs ----------
function wireDrop(dropId,inputId,onFile){const d=$(dropId),inp=$(inputId);
  inp.addEventListener('change',()=>{if(inp.files[0])onFile(inp.files[0])});
  d.addEventListener('dragover',e=>{e.preventDefault();d.classList.add('over')});d.addEventListener('dragleave',()=>d.classList.remove('over'));
  d.addEventListener('drop',e=>{e.preventDefault();d.classList.remove('over');const f=e.dataTransfer.files[0];if(f)onFile(f)})}
wireDrop('dropT','fileT',async f=>{setStatus('Reading '+f.name+' …');try{st.table=await readTable(f);st.tableName=f.name;
  const nS=st.table.kind==='sparse'?st.table.samples.length:Math.max(st.table.colIds.length,st.table.rowIds.length),nF=st.table.kind==='sparse'?st.table.feats.length:Math.min(st.table.colIds.length,st.table.rowIds.length);
  $('dropT').classList.add('done');$('dropT').querySelector('strong').textContent=f.name;$('dropT').querySelector('small').textContent='click to choose another file';$('infoT').innerHTML=`Read <b>${nF.toLocaleString()}</b> features and <b>${nS.toLocaleString()}</b> samples.`;stat('statT',true,`${nS.toLocaleString()} samples, ${nF.toLocaleString()} features`);setStatus('');ready()}
  catch(e){setStatus('Could not read the table: '+e.message,true)}});
wireDrop('dropM','fileM',async f=>{try{st.meta=await readMeta(f);st.metaName=f.name;$('dropM').classList.add('done');$('dropM').querySelector('strong').textContent=f.name;$('dropM').querySelector('small').textContent='click to choose another file';fillColumns();$('infoM').textContent=`${st.meta.rows.length} rows, ${st.meta.head.length-1} columns.`;stat('statM',true,`${st.meta.rows.length.toLocaleString()} rows, ${st.meta.head.length-1} columns`);ready()}
  catch(e){setStatus('Could not read the metadata: '+e.message,true)}});

// ---------- demo ----------
$('demoHero')&&$('demoHero').addEventListener('click',()=>$('demo').click());
$('demo').addEventListener('click',()=>{
  const lines=DEMO.table.split('\n'),feats=Array.from({length:DEMO.nFeat},(_,k)=>'ASV'+(k+1));
  const samples=lines.map(l=>{const [id,rest]=l.split('\t');let p=0;const idx=[],val=[];for(const t of rest.split(' ')){const [d,v]=t.split(':');p+=parseInt(d,36);idx.push(p);val.push(parseInt(v,36))}return {id,idx,val}});
  st.table={kind:'sparse',feats,samples};st.tableName='demo';stat('statT',true,`demo, ${samples.length} samples`);st.meta=parseMeta(DEMO.meta);st.metaName='demo metadata';
  $('dropT').classList.add('done');$('dropT').querySelector('strong').textContent='Demo table loaded';$('dropT').querySelector('small').textContent='or open your own table';$('infoT').innerHTML=`16 students, weekly gut samples (Student Microbiome Project, ${samples.length} samples). Two errors were planted, and the plate layout is simulated: can SwapSeq find them?`;
  $('dropM').classList.add('done');$('dropM').querySelector('strong').textContent='Demo metadata loaded';$('dropM').querySelector('small').textContent='or open your own sheet';fillColumns();ready();run();
});

// ---------- run in a worker ----------
let worker=null;
function getWorker(){if(worker)return worker;
  const src=$('coreSrc').textContent+`
  self.onmessage=e=>{const {groups,prior}=e.data;const out=[];
    for(const g of groups){const t0=Date.now();
      if(g.samples.length<4){out.push({name:g.name,skipped:'Fewer than 4 samples.'});continue}
      const data={feats:{length:g.nFeat},samples:g.samples.map(s=>({id:s.id,idx:Int32Array.from(s.idx),val:Float64Array.from(s.val)}))};
      const progress=m=>self.postMessage({type:'progress',msg:(groups.length>1?g.name+': ':'')+m});
      if(scorable(g.lab).filter(Boolean).length<4){out.push({name:g.name,skipped:'Fewer than 4 samples belong to people with two or more samples.'});continue}
      const M=distances(data,progress);const r=annotate(data,g.lab,analyse(data,g.lab,{prior,progress},M),{house:g.house,well:g.well,plate:g.plate});
      progress('Drawing the sample map');const map=g.samples.length<=3000?pcoa(M.J,g.samples.length,40):null;
      const pos=r.planted.pair.filter(Number.isFinite),obs=r.rows.map(x=>x.pair).filter(Number.isFinite);
      out.push({name:g.name,map:map?{x:Array.from(map.x),y:Array.from(map.y),share:map.share}:null,rows:r.rows,pairs:r.pairs,shifts:r.shifts||[],power:r.power,nSubjects:r.nSubjects,nScorable:r.nScorable,aitFeat:r.aitFeat,
        hist:{obs:r.rows.map(x=>x.single).filter(Number.isFinite),plant:r.planted.single,obsPair:obs,plantPair:pos},secs:(Date.now()-t0)/1000})}
    self.postMessage({type:'done',out})}`;
  worker=new Worker(URL.createObjectURL(new Blob([src],{type:'text/javascript'})));return worker}

$('run').addEventListener('click',run);
function run(){
  if(!st.table||!st.meta){const miss=!st.table&&!st.meta?'Open a feature table (step 1) and your sample metadata (step 2) first, or load the demo study.':!st.table?'Open a feature table in step 1 first.':'Open your sample metadata in step 2 first.';
    setStatus(miss,true);(!st.table?$('dropT'):$('dropM')).scrollIntoView({behavior:'smooth',block:'center'});return}const m=st.meta,kP=m.head.indexOf($('colP').value),gc=$('colG').value,kG=gc?m.head.indexOf(gc):-1,col=id=>$(id).value?m.head.indexOf($(id).value):-1,kH=col('colH'),kW=col('colW'),kL=col('colL');
  const all=samplesFrom(st.table,m);const groups=new Map();let skipped=0;
  for(const s of all.samples){const r=m.map.get(s.id);if(!r||!r[kP]||/^(na|nan|none|null|missing|not applicable)$/i.test(r[kP])){skipped++;continue}
    const g=kG>=0?(r[kG]||'(blank)'):'All samples';if(!groups.has(g))groups.set(g,{name:g,samples:[],lab:[],nFeat:all.feats,house:kH>=0?[]:null,well:kW>=0?[]:null,plate:kL>=0?[]:null});const G=groups.get(g);G.samples.push(s);G.lab.push(r[kP]);
    if(kH>=0)G.house.push(r[kH]||'');if(kW>=0)G.well.push(r[kW]||'');if(kL>=0)G.plate.push(r[kL]||'')}
  const list=[...groups.values()];const n=list.reduce((a,g)=>a+g.samples.length,0);
  if(!n){setStatus('No samples have a person ID in that column.',true);return}
  const big=list.find(g=>g.samples.length>4000);if(big){setStatus(`"${big.name}" has ${big.samples.length} samples. SwapSeq handles up to 4,000 per check in the browser; split it with "Check separately by".`,true);return}
  $('run').disabled=true;$('status').innerHTML='<span id="pmsg">Starting …</span><div class="progress"><i></i></div>';
  const w=getWorker();const prior=+$('prior').value;
  w.onmessage=e=>{if(e.data.type==='progress'){const p=$('pmsg');if(p)p.textContent=e.data.msg+' …'}else{ $('run').disabled=false;setStatus(`Checked ${n.toLocaleString()} samples${skipped?` (${skipped} without a person ID were left out)`:''}.`);render(e.data.out,prior)}};
  w.onerror=e=>{$('run').disabled=false;setStatus('Something went wrong: '+e.message,true)};
  st.lastCols={person:$('colP').value,group:gc,house:$('colH').value,well:$('colW').value,plate:$('colL').value};w.postMessage({groups:list,prior});
}

// ---------- results ----------
let lastOut=null;
function render(out,prior){lastOut=out;const box=$('groups');box.innerHTML='';$('results').hidden=false;
  $('resMeta').textContent=`${st.tableName}, expected error rate ${prior*100}%`;
  for(const g of out){const el=document.createElement('div');el.className='group';
    if(g.skipped){el.innerHTML=`<h3 style="margin-bottom:8px">${esc(g.name)}</h3><div class="verdict weak">Not checked: ${esc(g.skipped)}</div>`;box.appendChild(el);continue}
    const flagged=g.rows.filter(r=>r.level==='flag').sort((a,b)=>b.prob-a.prob),checks=g.rows.filter(r=>r.level==='check').sort((a,b)=>Math.max(b.broad,b.prob)-Math.max(a.broad,a.prob));
    const unscored=g.rows.filter(r=>!r.scorable).length,pw=g.power.swap,pw1=g.power.single,pwc=g.power.swapWithCheck,pw1c=g.power.singleWithCheck;
    const fPairs=g.pairs.filter(p=>p.level==='flag').length,cPairs=g.pairs.length-fPairs,nPeople=new Set(g.rows.filter(r=>r.scorable).map(r=>r.subject)).size,fewPeople=nPeople<8,lowPower=pwc<.25||fewPeople;
    let verdict;
    if(fewPeople)verdict=`<div class="verdict weak"><b>Too few people to calibrate.</b> Only ${nPeople} people have two or more samples. SwapSeq learns what an error looks like by swapping labels between people, and with fewer than 8 its power estimate is unreliable (in our benchmark, 5 mice gave an estimate of 36% when the real rate was 3%). Strong flags are still shown; the check list is hidden.</div>`;
    else if(lowPower)verdict=`<div class="verdict weak"><b>Low power in this dataset.</b> Even with the check list, only about ${fmtPct(pwc)} of planted swaps could be told apart from normal variation, so a clean result here does not prove the labels are right, and the check list below is mostly natural look-alikes. This is expected for infants, low-biomass sites and very shallow sequencing.</div>`;
    else if(flagged.length)verdict=`<div class="verdict warn"><b>${flagged.length} sample${flagged.length>1?'s':''} flagged${fPairs?`, including ${fPairs} likely swap${fPairs>1?'s':''}`:''}${checks.length?`, plus ${checks.length} to check`:''}.</b> Flags carry strong evidence. In this dataset SwapSeq flags about ${fmtPct(pw)} of swapped tubes, and ${fmtPct(pwc)} counting the check list.</div>`;
    else if(checks.length)verdict=`<div class="verdict weak"><b>No strong flags, ${checks.length} sample${checks.length>1?'s':''} to check.</b> These look somewhat like another person. Some will be natural look-alikes (people who live together, unusual samples), so compare them with lab records before changing anything.</div>`;
    else verdict=`<div class="verdict good"><b>No sample looks like it belongs to someone else.</b> In this dataset SwapSeq would have caught about ${fmtPct(pwc)} of swapped tubes and ${fmtPct(pw1c)} of one-off mislabels.</div>`;
    el.innerHTML=`${out.length>1?`<h3 style="margin-bottom:10px;font-size:19px">${esc(g.name)}</h3>`:''}
      <div class="tiles">
        <div class="panel tile"><div class="label">Samples checked</div><div class="v">${g.nScorable.toLocaleString()}</div><div class="d">${g.nSubjects} people${unscored?`, ${unscored} with a single sample not checked`:''}</div></div>
        <div class="panel tile"><div class="label">Flagged, and to check</div><div class="v"><span style="color:${flagged.length?'var(--swap)':'var(--ok)'}">${flagged.length}</span> <span style="color:var(--mis);font-size:.7em">+ ${checks.length}</span></div><div class="d">${fPairs} flagged swap pair${fPairs===1?'':'s'}${cPairs?`, ${cPairs} possible`:''}</div></div>
        <div class="panel tile"><div class="label">Catches swaps</div><div class="v">${fmtPct(pwc)}</div><div class="d">${fmtPct(pw)} as strong flags, estimated on errors planted in your data</div></div>
        <div class="panel tile"><div class="label">Catches mislabels</div><div class="v">${fmtPct(pw1c)}</div><div class="d">${fmtPct(pw1)} as strong flags, for one sample under the wrong person</div></div>
      </div>${verdict}${(g.shifts||[]).map(x=>`<div class="verdict warn"><b>Possible sample-sheet shift${x.plate?` on ${esc(x.plate)}`:''}:</b> ${x.n} tubes in a row, from position <span class="mono">${esc(x.from)}</span> to <span class="mono">${esc(x.to)}</span>, each look like the person written for the ${x.direction>0?'next':'previous'} tube. Check whether a row was skipped or added in the sample sheet.</div>`).join('')}`;
    const list=flagged.concat(lowPower?[]:checks);
    const tw=document.createElement('div');tw.className='tablewrap';
    tw.innerHTML=list.length?`<table><thead><tr><th>Sample</th><th>Labelled as</th><th>Verdict</th><th>Evidence</th><th>Looks like</th><th>Why</th></tr></thead><tbody>${list.map(r=>{
      const kind=r.kind,lvl=r.level,tag=r.cohab?'Housemate?':(kind==='swap'?'Swap':'Mismatch')+(lvl==='check'?'?':''),p=lvl==='flag'?r.prob:Math.max(r.broad,r.prob);
      let why=kind==='swap'?`Swapping labels with <span class="mono">${esc(r.partner)}</span> (labelled ${esc(r.partnerSubject)}) makes both samples fit their person.`
        :`Closest sample is <span class="mono">${esc(r.nearest)}</span> from ${esc(r.looksLike)}, not one of ${esc(r.subject)}'s own samples.`;
      if(r.rare&&r.rare.total)why+=`<span class="ev">Shares <b>${r.rare.withOther}</b> of its ${r.rare.total} rare ASVs with ${esc(r.rare.other)}'s other samples, and <b>${r.rare.withOwn}</b> with ${esc(r.subject)}'s.</span>`;
      if(r.shift)why+=`<span class="ev warnc">Part of a possible sample-sheet shift (positions ${esc(r.shift.from)} to ${esc(r.shift.to)}).</span>`;
      if(r.cohab)why+=`<span class="ev warnc">Same household (${esc(r.cohab)}) as ${esc(r.rare?r.rare.other:r.looksLike)}: people who live together share microbes${r.downgraded?', so this was moved to the check list':''}.</span>`;
      const pw=kind==='swap'?g.pairs.find(p=>(p.a===r.id||p.b===r.id)&&p.wells):null;
      if(pw){const w0=pw.wells,w=pw.a===r.id?w0:{...w0,a:w0.b,b:w0.a},rel={adjacent:'side by side',mirror:'mirror positions in one column (a flipped strip?)',diagonal:'diagonal neighbours',apart:'not neighbours'}[w.relation];why+=`<span class="ev">This tube is in well <span class="mono">${esc(w.a)}</span>, its partner in <span class="mono">${esc(w.b)}</span>${w.plate?` on ${esc(w.plate)}`:''}: ${rel}.</span>`}
      else if(r.well)why+=`<span class="ev">Well <span class="mono">${esc(r.well)}</span>.</span>`;
      return `<tr data-sid="${esc(r.id)}" tabindex="0" title="Show on the sample map"><td class="id">${esc(r.id)}</td><td class="mono">${esc(r.subject)}</td><td><span class="tag ${lvl==='flag'?kind:'review'}">${tag}</span></td><td><div class="prob"><span class="bar"><i style="width:${Math.round(p*100)}%;background:${lvl==='flag'?'var(--swap)':'var(--mis)'}"></i></span>${lvl==='flag'?'Strong':'Check'}</div></td><td class="mono">${esc(kind==='swap'?r.partnerSubject:r.looksLike)}</td><td class="why">${why}</td></tr>`}).join('')}</tbody></table>`
      :`<div class="empty">${lowPower&&checks.length?`${checks.length} weak look-alikes hidden because power is low here. They are in the copied results.`:'Nothing to show: every checked sample is closest to its own person.'}</div>`;
    el.appendChild(tw);
    const tools=document.createElement('div');tools.className='tools';
    tools.innerHTML=`<button class="btn ghost" data-act="copy">Copy QC audit (TSV)</button>`+(window.self===window.top?`<button class="btn ghost" data-act="dl">Download QC audit (.tsv)</button>`:'');
    tools.addEventListener('click',e=>{const a=e.target.dataset.act;if(!a)return;const tsv=toTSV(g);
      if(a==='copy'){navigator.clipboard.writeText(tsv).then(()=>{e.target.textContent='Copied'},()=>{fallbackCopy(tsv);e.target.textContent='Copied'})}
      else{const url=URL.createObjectURL(new Blob([tsv],{type:'text/tab-separated-values'}));const l=document.createElement('a');l.href=url;l.download=`swapseq_qc_${g.name.replace(/\W+/g,'_')}.tsv`;document.body.appendChild(l);l.click();l.remove();setTimeout(()=>URL.revokeObjectURL(url),2000)}});
    el.appendChild(tools);
    const ms=methodsText(g,prior);const mb=document.createElement('div');mb.className='panel methods';
    mb.innerHTML=`<div class="label">For your Methods section</div><p class="mtext">${esc(ms)}</p><button class="btn ghost" data-m="1">Copy text</button>`;
    mb.querySelector('button').addEventListener('click',e=>{navigator.clipboard.writeText(ms).then(()=>{e.target.textContent='Copied'},()=>{fallbackCopy(ms);e.target.textContent='Copied'})});
    el.appendChild(mb);
    if(g.map){const mp=sampleMap(g);el.insertBefore(mp,tools);
      tw.addEventListener('click',e=>{const tr=e.target.closest('tr[data-sid]');if(tr)mp.select(tr.dataset.sid,true)});
      tw.addEventListener('keydown',e=>{if(e.key==='Enter'){const tr=e.target.closest('tr[data-sid]');if(tr)mp.select(tr.dataset.sid,true)}})}
    el.appendChild(chart(g));
    box.appendChild(el)}
  $('results').scrollIntoView({behavior:'smooth',block:'start'});
}
function fallbackCopy(t){const ta=document.createElement('textarea');ta.value=t;document.body.appendChild(ta);ta.select();try{document.execCommand('copy')}catch(e){}ta.remove()}
const VERSION='SwapSeq method v3';
function methodsText(g,prior){const f=g.rows.filter(r=>r.level==='flag').length,c=g.rows.filter(r=>r.level==='check').length,fp=g.pairs.filter(p=>p.level==='flag').length;
  return `Sample identity was checked with ${VERSION} (soft Jaccard and Aitchison fingerprints, calibrated on errors planted in this dataset; expected mislabel rate ${prior*100}%). Of ${g.nScorable} samples from ${g.nSubjects} participants${g.name!=='All samples'?` (${g.name})`:''}, ${f} ${f===1?'was':'were'} flagged with strong evidence${fp?` (${fp} reciprocal swap${fp===1?'':'s'})`:''} and ${c} listed for checking against laboratory records. Estimated detection in this dataset: ${fmtPct(g.power.swapWithCheck)} of swapped tubes and ${fmtPct(g.power.singleWithCheck)} of single mislabels (${fmtPct(g.power.swap)} and ${fmtPct(g.power.single)} as strong flags).`}
function toTSV(g){const h=['sample','labelled_as','level','verdict','p_strict','p_broad','p_mislabel','p_swap','looks_like','nearest_sample','swap_partner','partner_labelled_as','rare_asvs','rare_shared_with_look_alike','rare_shared_with_own','well','same_household','mismatch_score','swap_score'];
  const f=v=>Number.isFinite(v)?(+v).toFixed(4):'',pairOf=id=>g.pairs.find(p=>p.a===id||p.b===id);
  const head=[`# ${VERSION} QC audit`,`# file: ${st.tableName}; group: ${g.name}; date: ${new Date().toISOString().slice(0,10)}`,`# columns used: ${Object.entries(st.lastCols||{}).filter(([k,v])=>v).map(([k,v])=>k+'='+v).join(', ')}`,`# ${methodsText(g,+$('prior').value)}`];
  return head.join('\n')+'\n'+h.join('\t')+'\n'+g.rows.map(r=>[r.id,r.subject,!r.scorable?'not_checked':r.level,r.level==='ok'?'ok':r.kind,f(r.prob),f(r.broad),f(r.pSingle),f(r.pSwap),r.looksLike??'',r.nearest??'',r.kind==='swap'?r.partner:'',r.kind==='swap'?r.partnerSubject:'',
    r.rare?r.rare.total:'',r.rare?r.rare.withOther:'',r.rare?r.rare.withOwn:'',r.well||'',r.cohab||'',f(r.single),f(r.pair)].join('\t')).join('\n')+'\n'}

// histogram: how this study's samples score against errors planted in it
function chart(g){const box=document.createElement('div');box.className='panel chart';
  const obs=g.hist.obsPair,pl=g.hist.plantPair;if(obs.length<5||pl.length<5){return box}
  const all=obs.concat(pl),lo=Math.min(...all),hi=Math.max(...all),B=40,w=(hi-lo)/B||1;
  const bin=x=>{const c=new Array(B).fill(0);for(const v of x)c[Math.min(B-1,Math.floor((v-lo)/w))]++;return c.map(v=>v/x.length)};
  const a=bin(obs),b=bin(pl),mx=Math.max(...a,...b);
  const W=900,H=230,L=44,R=12,T=10,Bt=36,iw=W-L-R,ih=H-T-Bt,xs=v=>L+(v-lo)/(hi-lo||1)*iw,ys=v=>T+ih-(v/mx)*ih;
  const flaggedIds=new Set(g.rows.filter(r=>r.flag).map(r=>r.id));
  let s=`<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Swap scores of your samples compared with planted swaps">`;
  for(let k=0;k<=4;k++){const y=T+ih*k/4;s+=`<line x1="${L}" x2="${W-R}" y1="${y}" y2="${y}" stroke="var(--grid)" stroke-width="1"/>`}
  a.forEach((v,k)=>{if(!v)return;const x=L+k*iw/B;s+=`<rect x="${x+1}" y="${ys(v)}" width="${iw/B-2}" height="${T+ih-ys(v)}" rx="2" fill="var(--obs)" opacity=".75"><title>${Math.round(v*obs.length)} of your samples, score ${(lo+k*w).toFixed(2)} to ${(lo+(k+1)*w).toFixed(2)}</title></rect>`});
  let d='';b.forEach((v,k)=>{const x0=L+k*iw/B,x1=x0+iw/B,y=ys(v);d+=(k?'L':'M')+x0+','+y+'L'+x1+','+y});s+=`<path d="${d}" fill="none" stroke="var(--plant)" stroke-width="2"/>`;
  if(lo<0&&hi>0)s+=`<line x1="${xs(0)}" x2="${xs(0)}" y1="${T}" y2="${T+ih}" stroke="var(--muted)" stroke-dasharray="4 4"/><text x="${xs(0)+6}" y="${T+14}" font-size="12" fill="var(--muted)" font-family="var(--sans)">swap would fit better →</text>`;
  g.rows.forEach(r=>{if(flaggedIds.has(r.id)&&Number.isFinite(r.pair))s+=`<circle cx="${xs(r.pair)}" cy="${T+ih-6}" r="5" fill="var(--swap)" stroke="var(--surface)" stroke-width="2"><title>${esc(r.id)}: flagged</title></circle>`});
  const ticks=5;for(let k=0;k<=ticks;k++){const v=lo+(hi-lo)*k/ticks;s+=`<text x="${xs(v)}" y="${H-14}" font-size="12" text-anchor="middle" fill="var(--muted)" font-family="var(--sans)">${v.toFixed(1)}</text>`}
  s+=`<line x1="${L}" x2="${W-R}" y1="${T+ih}" y2="${T+ih}" stroke="var(--line2)"/></svg>`;
  box.innerHTML=`<div class="label" style="margin-bottom:8px">Swap score: your samples vs. errors planted in your data</div><div class="legend"><span><i style="background:var(--obs)"></i>Your samples</span><span><i style="background:var(--plant)"></i>Planted swaps</span><span><i style="background:var(--swap);border-radius:50%"></i>Flagged</span></div>${s}<p class="hint">The further apart the two shapes are, the better SwapSeq can tell a swapped tube from normal variation in this study.</p>`;
  return box}

// ---------- validation table ----------



// ---------- sample map: PCoA of the fingerprint distances ----------
function sampleMap(g){
  const box=document.createElement('div');box.className='panel chart smap';
  const rows=g.rows,n=rows.length,X=g.map.x,Y=g.map.y;
  const W=760,H=470,L=24,R=24,T=20,Bm=36;
  const mnx=Math.min(...X),mxx=Math.max(...X),mny=Math.min(...Y),mxy=Math.max(...Y);
  const sx=v=>L+(v-mnx)/((mxx-mnx)||1)*(W-L-R),sy=v=>T+(1-(v-mny)/((mxy-mny)||1))*(H-T-Bm);
  const byPerson=new Map();rows.forEach((r,i)=>{if(!byPerson.has(r.subject))byPerson.set(r.subject,[]);byPerson.get(r.subject).push(i)});
  const idx=new Map(rows.map((r,i)=>[r.id,i]));
  const order=[...rows.keys()].sort((a,b)=>(rows[a].level!=='ok')-(rows[b].level!=='ok'));
  let dots='';for(const i of order){const r=rows[i],cls=r.level==='flag'?(r.kind==='swap'?'f-swap':'f-mis'):r.level==='check'?'f-check':'';
    dots+=`<circle class="pt ${cls}" data-i="${i}" cx="${sx(X[i]).toFixed(1)}" cy="${sy(Y[i]).toFixed(1)}" r="${r.level==='ok'?4:6.5}"/>`}
  const pct=v=>Math.round(v*100);
  box.innerHTML=`<div class="map-h"><h3>Sample map</h3><p class="hint">Each dot is a sample. Samples with similar microbial fingerprints sit close together, so one person's samples usually form a cluster. Point at a dot to light up everyone else from that person; click a flagged one to see where it really belongs.</p></div>
    <div class="legend"><span><i style="background:var(--obs)"></i>Sample</span><span><i style="background:var(--swap);border-radius:50%"></i>Flagged swap</span><span><i style="background:var(--mis);border-radius:50%"></i>Flagged or to check</span><span><i style="background:var(--accent);border-radius:50%"></i>Same person as the one you point at</span><span><i style="background:none;border:2px solid var(--ok);border-radius:50%"></i>Person it looks like</span></div>
    <div class="map-wrap"><svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Map of samples by microbial fingerprint">
      <line x1="${L}" x2="${W-R}" y1="${H-Bm+6}" y2="${H-Bm+6}" stroke="var(--line)"/>
      <text x="${L}" y="${H-10}" font-size="12" fill="var(--faint)" font-family="var(--sans)">Horizontal axis ${pct(g.map.share[0])}% of variation, vertical ${pct(g.map.share[1])}%; PCoA of soft-Jaccard distances</text>
      <g class="links"></g><g class="dots">${dots}</g></svg><div class="tip" hidden></div></div>
    <p class="map-sel" aria-live="polite"></p>`;
  const svg=box.querySelector('svg'),tip=box.querySelector('.tip'),sel=box.querySelector('.map-sel'),links=box.querySelector('.links'),pts=[...box.querySelectorAll('.pt')];
  const ptOf=i=>pts.find(p=>+p.dataset.i===i);
  const clear=()=>pts.forEach(p=>p.classList.remove('hl','own','other','dim','sel'));
  let pinned=null;
  const light=(i,pin)=>{clear();links.innerHTML='';const r=rows[i];
    pts.forEach(p=>{const j=+p.dataset.i;if(rows[j].subject===r.subject)p.classList.add(j===i?'sel':'hl');else p.classList.add('dim')});
    if(pin&&r.level!=='ok'){const other=r.kind==='swap'?r.partnerSubject:r.looksLike;
      (byPerson.get(other)||[]).forEach(j=>{const p=ptOf(j);p.classList.remove('dim');p.classList.add('other')});
      const near=idx.get(r.kind==='swap'?r.partner:r.nearest);
      if(near!=null){links.innerHTML=`<line x1="${sx(X[i])}" y1="${sy(Y[i])}" x2="${sx(X[near])}" y2="${sy(Y[near])}" stroke="var(--ok)" stroke-width="2" stroke-dasharray="4 4"/>`}
      sel.innerHTML=`<b>${esc(r.id)}</b> is labelled <b>${esc(r.subject)}</b> (blue), but sits with <b>${esc(other)}</b> (green ring)${r.kind==='swap'?`. Its swap partner is <span class="mono">${esc(r.partner)}</span>.`:'.'} <button class="linkbtn" type="button">Clear</button>`;
      sel.querySelector('button').onclick=()=>{pinned=null;clear();links.innerHTML='';sel.textContent=''}}
  };
  svg.addEventListener('pointermove',e=>{const p=e.target.closest('.pt');if(!p){tip.hidden=true;if(pinned==null){clear()}else light(pinned,true);return}
    const i=+p.dataset.i,r=rows[i];if(pinned==null)light(i,false);
    const status=r.level==='flag'?(r.kind==='swap'?`Flagged swap with ${r.partner}`:`Flagged: looks like ${r.looksLike}`):r.level==='check'?`To check: looks like ${r.looksLike}`:'Fits its person';
    tip.innerHTML=`<b>${esc(r.id)}</b><br>Labelled ${esc(r.subject)}<br>${esc(status)}`;tip.hidden=false;
    const bb=box.querySelector('.map-wrap').getBoundingClientRect(),x=e.clientX-bb.left,y=e.clientY-bb.top;
    tip.style.left=Math.min(bb.width-190,x+14)+'px';tip.style.top=Math.max(0,y-58)+'px'});
  svg.addEventListener('pointerleave',()=>{tip.hidden=true;if(pinned==null)clear();else light(pinned,true)});
  svg.addEventListener('click',e=>{const p=e.target.closest('.pt');if(!p)return;const i=+p.dataset.i;pinned=rows[i].level!=='ok'?i:null;light(i,pinned!=null)});
  box.select=(id,scroll)=>{const i=idx.get(id);if(i==null)return;pinned=i;light(i,true);if(scroll)box.scrollIntoView({behavior:'smooth',block:'center'})};
  return box;
}

// ---------- validation charts (public benchmark) ----------
(function(){if(typeof VAL==='undefined')return;
  const ORDER=['SMP gut','SMP forehead','SMP palm','IMP adults (16S)','SMP tongue','SMP gut (2 per person)','SMP forehead (2 per person)','Diet study (species)','Diet study (strain-level)','Diet study (strain-level) (2 per person)','ECAM infants','MDSINE2 mice (healthy donor)','MDSINE2 mice (uc donor)'];
  const NAME={'IMP adults (16S)':'IMP gut, monthly','SMP gut (2 per person)':'SMP gut, 2 per person','SMP forehead (2 per person)':'SMP forehead, 2 per person','Diet study (species)':'Diet study, species','Diet study (strain-level)':'Diet study, strain level','Diet study (strain-level) (2 per person)':'Diet study, 2 per person','MDSINE2 mice (healthy donor)':'Mice, healthy donor','MDSINE2 mice (uc donor)':'Mice, UC donor'};
  const grp=d=>/2 per person/.test(d)?(/Diet/.test(d)?'Shotgun, adults sampled daily':'Before/after, 2 samples per person'):/^(SMP|IMP)/.test(d)?'16S, adults':/^Diet/.test(d)?'Shotgun, adults sampled daily':'Little individual signal (controls)';
  const by=new Map(VAL.ours.map((o,k)=>[o.dataset,[o,VAL.bc[k]]]));
  const rows=ORDER.filter(d=>by.has(d)).map(d=>{const [o,b]=by.get(d);return {d,name:NAME[d]||d,group:grp(d),n:o.samples,people:o.people,
    swaps:[o.swappedSamplesCaught,o.withReview.swaps,b.swappedSamplesCaught],mis:[o.oneWayCaught,o.withReview.oneWay,b.oneWayCaught],fa:[o.falseFlagsPer1000,o.withReview.falsePer1000,b.falseFlagsPer1000]}});
  const is16=r=>r.group==='16S, adults'||r.group==='Before/after, 2 samples per person';
  const rng=(k,i)=>{const v=rows.filter(is16).map(r=>r[k][i]).filter(x=>x!=null);return [Math.min(...v),Math.max(...v)]};
  const P=v=>Math.round(v*100)+'%';
  const M={
    swaps:{title:()=>{const a=rng('swaps',0),b=rng('swaps',2);return `On 16S data SwapSeq flags ${Math.round(a[0]*100)}–${P(a[1])} of swapped samples, Bray-Curtis ${Math.round(b[0]*100)}–${P(b[1])}`},max:1,fmt:P,unit:'of swapped samples caught',ticks:[0,.25,.5,.75,1]},
    mis:{title:()=>'Single mislabels are harder: SwapSeq flags up to 83%, and the check list adds more',max:1,fmt:P,unit:'of one-way mislabels caught',ticks:[0,.25,.5,.75,1]},
    fa:{title:()=>`Strong flags give at most ${Math.max(...rows.map(r=>r.fa[0])).toFixed(1)} false alarms per 1,000 samples in every dataset`,max:16,fmt:v=>Number.isInteger(v)?String(v):v.toFixed(1),unit:'false alarms per 1,000 untouched samples',ticks:[0,4,8,12,16]}};
  const host=document.getElementById('valChart');if(!host)return;let cur='swaps';
  const draw=()=>{const m=M[cur],x0=270,x1=736,y0=60,rh=26,gh=34,X=v=>x0+Math.min(1,v/m.max)*(x1-x0);
    let y=y0,g=null,parts=[];
    for(const r of rows){if(r.group!==g){if(g!==null)y+=gh;g=r.group;parts.push(`<text x="0" y="${y-22}" font-weight="700" fill="var(--ink)">${esc(g)}</text>`)}
      const [a,b,c]=r[cur];const lbl=`<text x="14" y="${y+4}" fill="var(--ink)">${esc(r.name)}</text>`;
      if(a==null){parts.push(lbl+`<text x="${x0}" y="${y+4}" fill="var(--faint)">not tested (two samples per person)</text>`);y+=rh;continue}
      const lo=Math.min(a,b,c),hi=Math.max(a,b,c);
      parts.push(`<g class="vrow" data-name="${esc(r.name)}" data-a="${m.fmt(a)}" data-b="${m.fmt(b)}" data-c="${m.fmt(c)}" data-n="${r.n}" data-p="${r.people}">${lbl}<rect x="0" y="${y-12}" width="760" height="24" fill="transparent"/><line x1="${X(lo)}" x2="${X(hi)}" y1="${y}" y2="${y}" stroke="var(--line2)" stroke-width="2.5" stroke-linecap="round"/>`+
        `<circle cx="${X(c)}" cy="${y}" r="5.5" fill="var(--obs)"/><circle cx="${X(b)}" cy="${y}" r="6.5" fill="var(--surface)" stroke="var(--accent)" stroke-width="2.2"/><circle cx="${X(a)}" cy="${y}" r="6.5" fill="var(--accent)"/></g>`);y+=rh}
    const H=y+30;
    const ax=m.ticks.map(t=>`<line x1="${X(t)}" x2="${X(t)}" y1="${y0-26}" y2="${y-6}" stroke="var(--grid)"/><text x="${X(t)}" y="${y+14}" text-anchor="middle" fill="var(--faint)">${m.fmt(t)}</text>`).join('');
    const lg=cur==='fa'?['Bray-Curtis strong flags','SwapSeq strong flags','SwapSeq flags + check list']:['Bray-Curtis, same pipeline','SwapSeq strong flag','SwapSeq flag or check list'];
    host.querySelector('.vtitle').textContent=m.title();host.querySelector('.vunit').textContent=`Shown: ${m.unit}. Point at a row for its numbers.`;
    host.querySelector('.vlegend').innerHTML=`<span><i style="background:var(--obs);border-radius:50%"></i>${lg[0]}</span><span><i style="background:var(--accent);border-radius:50%"></i>${lg[1]}</span><span><i style="border:2px solid var(--accent);border-radius:50%"></i>${lg[2]}</span>`;
    host.querySelector('.vsvg').innerHTML=`<svg viewBox="0 0 760 ${H}" role="img" aria-label="${esc(m.title())}" font-family="var(--sans)" font-size="13">${ax}${parts.join('')}</svg>`};
  host.innerHTML=`<div class="vtabs" role="group" aria-label="What to show"><button type="button" data-m="swaps" aria-pressed="true">Swaps caught</button><button type="button" data-m="mis" aria-pressed="false">Mislabels caught</button><button type="button" data-m="fa" aria-pressed="false">False alarms</button></div>
    <h3 class="vtitle"></h3><p class="hint vunit"></p><div class="legend vlegend"></div><div class="map-wrap"><div class="vsvg"></div><div class="tip" hidden></div></div>`;
  host.querySelector('.vtabs').addEventListener('click',e=>{const b=e.target.closest('button');if(!b)return;cur=b.dataset.m;host.querySelectorAll('.vtabs button').forEach(x=>x.setAttribute('aria-pressed',x===b));draw()});
  const tip=host.querySelector('.tip'),wrap=host.querySelector('.map-wrap');
  wrap.addEventListener('pointermove',e=>{const r=e.target.closest('.vrow');if(!r){tip.hidden=true;return}
    const fa=cur==='fa';tip.innerHTML=`<b>${r.dataset.name}</b><br>${r.dataset.n} samples, ${r.dataset.p} people<br>SwapSeq: ${r.dataset.a} strong${fa?'':', '+r.dataset.b+' with check list'}${fa?', '+r.dataset.b+' with check list':''}<br>Bray-Curtis: ${r.dataset.c}`;tip.hidden=false;
    const bb=wrap.getBoundingClientRect();tip.style.left=Math.min(bb.width-220,e.clientX-bb.left+14)+'px';tip.style.top=Math.max(0,e.clientY-bb.top-70)+'px'});
  wrap.addEventListener('pointerleave',()=>tip.hidden=true);
  draw();
  // plate-pattern errors on a real plating order
  const ph=document.getElementById('plateChart');if(!ph||!VAL.plate)return;
  const pr=VAL.plate,W=760,x0=250,x1=560,X=v=>x0+v*(x1-x0);let y=70,out='';
  const items=[['Neighbouring tubes swapped',pr.neighbour],['Sheet slips by 3 tubes',pr.L3],['Sheet slips by 5 tubes',pr.L5],['Sheet slips by 8 tubes',pr.L8]];
  for(const [lab,r] of items){const a=r.samplesFlagged,b=r.withCheck,c=r.shiftNamed;
    out+=`<g class="vrow" data-name="${lab}" data-a="${P(a)}" data-b="${P(b)}" data-c="${c==null?'–':P(c)}" data-f="${r.falseFlagsPer1000.toFixed(1)}"><rect x="0" y="${y-14}" width="760" height="40" fill="transparent"/><text x="0" y="${y+4}" fill="var(--ink)">${lab}</text>
      <rect x="${x0}" y="${y-9}" width="${X(b)-x0}" height="18" rx="4" fill="var(--accent)" fill-opacity=".25"/><rect x="${x0}" y="${y-9}" width="${X(a)-x0}" height="18" rx="4" fill="var(--accent)"/>
      <text x="${X(b)+10}" y="${y+4}" fill="var(--ink)" font-weight="600">${P(a)}</text><text x="${X(b)+46}" y="${y+4}" fill="var(--muted)">${P(b)} with check list</text>${c!=null?`<text x="0" y="${y+22}" fill="var(--muted)" font-size="12.5">shift named as a shift in ${P(c)} of runs</text>`:''}</g>`;y+=c!=null?56:46}
  ph.innerHTML=`<h3>Lab-style plate errors are caught too</h3><p class="hint">Errors planted on the real plating order of the diet study (526 tubes, 6 plates), 30 repeats each. No false shift report in any repeat.</p>
    <div class="legend"><span><i style="background:var(--accent)"></i>Strong flag</span><span><i style="background:var(--accent);opacity:.3"></i>Added by the check list</span></div>
    <div class="map-wrap"><svg viewBox="0 0 ${W} ${y}" role="img" aria-label="Plate-pattern errors caught" font-family="var(--sans)" font-size="13">${[0,.5,1].map(t=>`<line x1="${X(t)}" x2="${X(t)}" y1="44" y2="${y-20}" stroke="var(--grid)"/><text x="${X(t)}" y="40" text-anchor="middle" fill="var(--faint)">${P(t)}</text>`).join('')}${out}</svg><div class="tip" hidden></div></div>`;
  const t2=ph.querySelector('.tip'),w2=ph.querySelector('.map-wrap');
  w2.addEventListener('pointermove',e=>{const r=e.target.closest('.vrow');if(!r){t2.hidden=true;return}t2.innerHTML=`<b>${r.dataset.name}</b><br>${r.dataset.a} of planted errors flagged<br>${r.dataset.b} with the check list<br>Shift named: ${r.dataset.c}<br>False flags: ${r.dataset.f} per 1,000`;t2.hidden=false;
    const bb=w2.getBoundingClientRect();t2.style.left=Math.min(bb.width-220,e.clientX-bb.left+14)+'px';t2.style.top=Math.max(0,e.clientY-bb.top-80)+'px'});
  w2.addEventListener('pointerleave',()=>t2.hidden=true);
})();

