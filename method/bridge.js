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
