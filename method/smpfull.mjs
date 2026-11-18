// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mankirat Singh
// SMP full table (all samples, not only the time-series subset): does SwapSeq flag the two gut samples the SMP authors flagged?
import {analyse,annotate} from './core.mjs';import {biomJSON,readTSV} from './load.mjs';import fs from 'fs';import zlib from 'zlib';
const B=(process.env.OUT||'.')+'/';
const all=biomJSON(zlib.gunzipSync(fs.readFileSync(B+'data/smp/otu-tables/otu_table_even10000.biom.gz')).toString());
const m=new Map(readTSV(B+'data/smp/analysis-results/issue-13/map.txt').map(r=>[r['#SampleID'],r]));
for(const site of ['gut','tongue','palm','forehead']){
  let S=all.samples.filter(s=>m.has(s.id)&&m.get(s.id).BodySite===site&&!/^control/i.test(m.get(s.id).PersonalID));
  const c=new Map();S.forEach(s=>{const p=m.get(s.id).PersonalID;c.set(p,(c.get(p)||0)+1)});
  const data={feats:all.feats,samples:S},lab=S.map(s=>m.get(s.id).PersonalID);
  const r=annotate(data,lab,analyse(data,lab));const authors=S.filter(s=>m.get(s.id)['mislabeled_probability_above_0.90']==='TRUE').map(s=>s.id);
  console.log(`SMP ${site} (full table): ${S.length} samples, ${c.size} people; strong flags ${r.rows.filter(x=>x.level==='flag').length}, checks ${r.rows.filter(x=>x.level==='check').length}; authors flagged: ${authors.join(',')||'none'}`);
  for(const id of authors){const x=r.rows.find(y=>y.id===id);console.log('   author-flagged',id,x.subject,'nOwn',c.get(x.subject),'-> SwapSeq',x.level,x.kind,'looks like',x.looksLike,'strict',x.prob.toFixed(3),'broad',x.broad.toFixed(3),'single',x.single?.toFixed(3),x.rare?`rare ${x.rare.withOther}/${x.rare.total} own ${x.rare.withOwn}`:'')}
}
