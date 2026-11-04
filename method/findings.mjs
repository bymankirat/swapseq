// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mankirat Singh
import {analyse,annotate} from './core.mjs';import {smp,ecam,dietstudy,mdsine} from './load.mjs';
for(const ds of [smp('gut'),smp('tongue'),smp('palm'),smp('forehead'),ecam(),dietstudy('burst'),mdsine('healthy'),mdsine('uc')]){
  const r=annotate(ds.data,ds.lab,analyse(ds.data,ds.lab));const f=r.rows.filter(x=>x.level==='flag'),c=r.rows.filter(x=>x.level==='check');
  console.log(`${ds.name}: ${f.length} flags, ${c.length} checks, pairs ${r.pairs.length} (flag-level ${r.pairs.filter(p=>p.level==='flag').length})`);
  for(const x of f)console.log('   FLAG',x.kind,x.id,x.subject,'->',x.kind==='swap'?x.partnerSubject+' w/ '+x.partner:x.looksLike,'t',ds.time[x.i],x.rare?`rare ${x.rare.withOther}/${x.rare.total} own ${x.rare.withOwn}`:'')}
