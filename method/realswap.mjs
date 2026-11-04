// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mankirat Singh
// Blind test on a real, author-corrected tube swap (diet study): restore the original labels and run SwapSeq.
import {analyse,annotate,distances} from './core.mjs';import {dietstudy} from './load.mjs';
for(const lvl of ['burst','species']){const ds=dietstudy(lvl);const lab=ds.lab.slice(),ids=ds.data.samples.map(s=>s.id);
  const a=ids.indexOf('MCT.f.0086'),b=ids.indexOf('MCT.f.0417');[lab[a],lab[b]]=[lab[b],lab[a]]; // original pre-correction labels
  const M=distances(ds.data);
  for(const run of ['original labels (swap present)','corrected labels']){const L=run.startsWith('orig')?lab:ds.lab;
    const r=annotate(ds.data,L,analyse(ds.data,L,{},{...M,masked:false}));const fl=r.rows.filter(x=>x.level!=='ok');
    console.log(`${ds.name} | ${run}: ${fl.length} non-ok; pairs ${JSON.stringify(r.pairs.map(p=>[p.a,p.b,p.level]))}`);
    for(const x of fl)console.log('   ',x.level,x.kind,x.id,'labelled',x.subject,'->',x.kind==='swap'?x.partnerSubject:x.looksLike,'strict',x.prob.toFixed(3),'broad',x.broad.toFixed(3),x.rare?`rare ${x.rare.withOther}/${x.rare.total} vs own ${x.rare.withOwn}`:'')}}
