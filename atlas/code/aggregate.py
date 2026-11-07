import json,glob,os
R=[json.load(open(f)) for f in sorted(glob.glob('results/*.json'))]
GROUP={'smpGut':'Adult gut','diet':'Adult gut','imp':'Adult gut','david':'Adult gut','hagan':'Adult gut','mars':'Adult gut','poyet':'Adult gut','ibdmdb':'Adult gut','palleja':'Adult gut',
 'taur':'Hospital patients','smpTongue':'Oral','smpPalm':'Skin','smpForehead':'Skin','digiulio':'Vaginal',
 'ecam':'Infant gut','karelia':'Infant gut','shao':'Infant gut','heInfants':'Infant gut','kostic':'Infant gut','wandro':'Infant gut',
 'mdsH':'Mouse gut','mdsU':'Mouse gut','muVRE':'Mouse gut','storage':'Same-stool aliquots','dietOriginal':'Adult gut'}
rows=[]
for r in R:
    k=r['key'];g=GROUP.get(k,'Oral' if k.startswith('ploeg') else '?')
    m=r.get('measured') or {}
    c=m.get('swapCaught')
    grade=None if c is None else ('High' if c>=.8 else 'Medium' if c>=.5 else 'Low' if c>=.2 else 'Very low')
    P=json.load(open('private/'+k+'.json'));sup=sum(1 for q in P if q['level']=='flag' and q['rare'] and q['rare']['withOther']>=3 and q['rare']['withOther']>=2*max(1,q['rare']['withOwn']))
    rows.append(dict(key=k,supported=sup,name=r['name'],ref=r.get('ref'),type=r.get('type'),group=g,design=r.get('design'),validation=r['validation'],control=r['control'],
      samples=r['samples'],dropped=r['droppedLowDepth'],subjects=r['subjects'],scorable=r['scorable'],perSubject=r['medianSamplesPerSubject'],
      reid=round(r['reidentified']*100,1) if r['reidentified'] is not None else None,
      caught=None if c is None else round(c*100),review=None if not m else round(m['swapCaughtWithReview']*100),
      pairNamed=None if not m else round(m['swapPairNamed']*100),relabel=None if not m or m['relabelCaught'] is None else round(m['relabelCaught']*100),
      fp=None if not m else round(m['falseFlagsPer1000'],2),fpCheck=None if not m else round(m['falseChecksPer1000'],1),
      reps=m.get('reps'),testsPerRep=(m.get('swapsPerRep',0)*2+m.get('relabelsPerRep',0)) if m else None,
      flags=r['strongFlags'],checks=r['checkList'],pairs=r['strongSwapPairs'],allPairs=r['allSwapPairs'],
      flagsPer1000=round(r['flagsPer1000'],1) if r['flagsPer1000'] is not None else None,grade=grade,
      selfPowerSwap=round(r.get('selfEstimatedPowerOnOriginal',r.get('power'))['swap']*100) if r.get('selfEstimatedPowerOnOriginal',r.get('power'))['swap']==r.get('selfEstimatedPowerOnOriginal',r.get('power'))['swap'] else None))
rows.sort(key=lambda x:(-(x['caught'] if x['caught'] is not None else -1),x['name']))
main=[x for x in rows if not x['control']]
studies=set()
for x in main: studies.add(x['ref'].split(',')[0] if x['ref'] else x['name'])
S=dict(datasets=len(main),studies=len(studies),samples=sum(x['samples'] for x in main),subjects=sum(x['subjects'] for x in main),
  grades={g:sum(1 for x in main if x['grade']==g) for g in ['High','Medium','Low','Very low']},
  flags=sum(x['flags'] for x in main),checks=sum(x['checks'] for x in main),pairs=sum(x['pairs'] for x in main),
  datasetsWithFlags=sum(1 for x in main if x['flags']>0),supported=sum(x['supported'] for x in main),datasetsWithSupported=sum(1 for x in main if x['supported']>0))
json.dump({'summary':S,'rows':rows},open('atlas_data.json','w'),indent=1)
print(json.dumps(S,indent=1))
for x in rows: print(f"{x['key']:32}{x['group']:20}{x['samples']:6}{x['subjects']:5} reid {x['reid']!s:6} caught {x['caught']!s:5} rev {x['review']!s:5} rel {x['relabel']!s:5} fp {x['fp']!s:6} flags {x['flags']:3} sup {x['supported']} pairs {x['pairs']} {x['grade']} selfpow {x['selfPowerSwap']}")
