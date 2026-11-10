# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mankirat Singh
# Builds the Atlas table from the per-dataset results of run.mjs and writes it into the page.
#   python3 atlas/code/build_atlas.py RUN_DIR      (RUN_DIR holds results/*.json and private/*.json from run.mjs)
# Datasets whose source files could not be fetched again keep their row from data/atlas.json (see CARRY).
import json,glob,os,sys,csv,io,re
RUN=sys.argv[1];HERE=os.path.dirname(os.path.abspath(__file__));AT=os.path.dirname(HERE)
GROUP={'smpGut':'Adult gut','diet':'Adult gut','imp':'Adult gut','david':'Adult gut','hagan':'Adult gut','mars':'Adult gut','poyet':'Adult gut','ibdmdb':'Adult gut','palleja':'Adult gut',
 'taur':'Hospital patients','smpTongue':'Oral','smpPalm':'Skin and nose','smpForehead':'Skin and nose','digiulio':'Vaginal',
 'ecam':'Infant gut','karelia':'Infant gut','shao':'Infant gut','heInfants':'Infant gut','kostic':'Infant gut','wandro':'Infant gut',
 'yeBrooks':'Infant gut','yeT2d':'Adult gut','yeHall':'Gut, other cohorts','yePRJEB38984':'Gut, other cohorts','yePRJEB70966':'Gut, other cohorts','yeFMT':'Gut, other cohorts',
 'mdsH':'Mouse gut','mdsU':'Mouse gut','muVRE':'Mouse gut','storage':'Same-stool aliquots','dietOriginal':'Adult gut',
 'hmpStool':'Adult gut','cafe':'Adult gut','hmpTongue':'Oral','hmpPlaque':'Oral','hmpSkin':'Skin and nose','hmpNostril':'Skin and nose'}
SHORT={'hmpStool':'HMP1 stool','hmpTongue':'HMP1 tongue','hmpPlaque':'HMP1 dental plaque','hmpSkin':'HMP1 skin behind the ear','hmpNostril':'HMP1 nostril','cafe':'Controlled feeding (Wu)'}
SKIP={'turnbaugh'}   # two samples per person: the Atlas test needs people with three or more
old={r['key']:r for r in json.load(open(AT+'/data/atlas.json'))}
rows=[]
for f in sorted(glob.glob(RUN+'/results/*.json')):
    r=json.load(open(f));k=r['key']
    if k in SKIP:continue
    m=r.get('measured') or {};c=m.get('swapCaught')
    grade=None if c is None else ('High' if c>=.8 else 'Medium' if c>=.5 else 'Low' if c>=.2 else 'Very low')
    P=json.load(open(RUN+'/private/'+k+'.json'))
    sup=sum(1 for q in P if q['level']=='flag' and q['rare'] and q['rare']['withOther']>=3 and q['rare']['withOther']>=2*max(1,q['rare']['withOwn']))
    sp=r.get('selfEstimatedPowerOnOriginal',r.get('power'))['swap'];rp=r.get('reportedPower')
    rows.append(dict(key=k,supported=sup,name=r['name'],ref=r.get('ref'),type=r.get('type'),group=GROUP.get(k,'Oral' if k.startswith('ploeg') else '?'),design=r.get('design'),
      validation=r['validation'],control=r['control'],samples=r['samples'],dropped=r['droppedLowDepth'],subjects=r['subjects'],scorable=r['scorable'],perSubject=r['medianSamplesPerSubject'],
      reid=round(r['reidentified']*100,1) if r['reidentified'] is not None else None,caught=None if c is None else round(c*100),
      review=None if not m else round(m['swapCaughtWithReview']*100),pairNamed=None if not m else round(m['swapPairNamed']*100),
      relabel=None if not m or m['relabelCaught'] is None else round(m['relabelCaught']*100),fp=None if not m else round(m['falseFlagsPer1000'],2),
      fpCheck=None if not m else round(m['falseChecksPer1000'],1),reps=m.get('reps'),testsPerRep=(m.get('swapsPerRep',0)*2+m.get('relabelsPerRep',0)) if m else None,
      flags=r['strongFlags'],checks=r['checkList'],pairs=r['strongSwapPairs'],allPairs=r['allSwapPairs'],
      flagsPer1000=round(r['flagsPer1000'],1) if r['flagsPer1000'] is not None else None,grade=grade,selfPowerSwap=round(sp*100) if sp==sp else None,
      short=old[k]['short'] if k in old else SHORT[k],powerV111=None if not rp else round(100*rp['swap']),possibleDuplicates=r['possibleDuplicates'],engine='November update 2'))
have={r['key'] for r in rows};CARRY=[k for k in old if k not in have]
for k in CARRY:rows.append({**old[k],'group':GROUP.get(k,old[k]['group']),'engine':'November update (source files not available for this run)'})
rows.sort(key=lambda x:(-(x['caught'] if x['caught'] is not None else -1),x['name']))
json.dump(rows,open(AT+'/data/atlas.json','w'),indent=1)
keys=list(rows[0].keys());buf=io.StringIO();w=csv.DictWriter(buf,fieldnames=keys);w.writeheader();w.writerows(rows);open(AT+'/data/atlas.csv','w').write(buf.getvalue())
page=open(AT+'/index.html').read();a=page.index('<script>const ROWS=[');b=page.index('\n',a)
line=page[a:b];tail=line[line.rindex('];'):]   # keep whatever follows the array on that line
page=page[:a]+'<script>const ROWS='+json.dumps(rows,separators=(',',':'),ensure_ascii=False)[:-1]+tail+page[b:]
open(AT+'/index.html','w').write(page)
main=[x for x in rows if not x['control']]
print('datasets',len(main),'samples',sum(x['samples'] for x in main),'people',sum(x['subjects'] for x in main),'studies',len({(x['ref'] or x['name']).split(',')[0].split(';')[0] for x in main}),
 'grades',{g:sum(1 for x in main if x['grade']==g) for g in ['High','Medium','Low','Very low']},'>=50:',sum(1 for x in main if (x['caught'] or 0)>=50),'carried over:',CARRY)
