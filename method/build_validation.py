# Builds the public Methods & validation page from the same data and chart code as the tool.
import json,re
B='/tmp/claude-0/-home-claude/d674b933-5c36-5fd7-9c63-4fe2a591bf95/scratchpad/swap/'
tool=open(B+'tool.html').read();app=open(B+'app.js').read();t=open(B+'validation_src.html').read()
i=tool.index('<style>');css=tool[i+7:tool.index('</style>',i)]
a=app.index('// ---------- validation charts');b=app.index('\n})();',a)+6;charts=app[a:b]
L=lambda f:json.load(open(B+f))
val={'ours':L('bench_results_prior0.01_v3.json')+L('bench_results_prior0.01_v3new.json')+L('bench_results_prior0.01_imp.json'),
     'bc':L('bench_results_prior0.01_braycurtis_v3.json')+L('bench_results_prior0.01_braycurtis_v3new.json')+L('bench_results_prior0.01_braycurtis_imp.json')}
assert [o['dataset'] for o in val['ours']]==[o['dataset'] for o in val['bc']]
plate={k:L(f'plate_results_burst_{f}.json')['res'][n] for k,f,n in [('neighbour','L5','neighbour swap'),('L3','L3','shift run'),('L5','L5','shift run'),('L8','L8','shift run')]}
plate['L3']['shiftNamed']=plate['L3'].get('shiftNamed',0);plate['neighbour']['shiftNamed']=None
pc=lambda v:'–' if v is None else '%d%%'%round(v*100)
rws=''
for o,bc in zip(val['ours'],val['bc']):
  w=o['withReview'];rws+=f"<tr><td>{o['dataset']}</td><td class='num'>{o['samples']}</td><td class='num'>{o['people']}</td><td class='num best'>{pc(o['swappedSamplesCaught'])}</td><td class='num'>{pc(w['swaps'])}</td><td class='num'>{pc(bc['swappedSamplesCaught'])}</td><td class='num'>{pc(o['swapPairsNamedExactly'])}</td><td class='num'>{pc(o['oneWayCaught'])}</td><td class='num'>{pc(w['oneWay'])}</td><td class='num'>{pc(bc['oneWayCaught'])}</td><td class='num'>{o['falseFlagsPer1000']:.1f}</td><td class='num'>{w['falsePer1000']:.1f}</td><td class='num'>{pc(o['selfEstimatedPower']['swap'])}</td></tr>"
th="<thead><tr><th>Dataset</th><th class='num'>Samples</th><th class='num'>People</th><th class='num'>Swaps flagged</th><th class='num'>+ check list</th><th class='num'>Bray-Curtis</th><th class='num'>Pairs named</th><th class='num'>Mislabels flagged</th><th class='num'>+ check list</th><th class='num'>Bray-Curtis</th><th class='num'>False flags / 1000</th><th class='num'>False checks / 1000</th><th class='num'>Reported power</th></tr></thead>"
t=t.replace('<table class="val" id="valTable"></table>','<table class="val" id="valTable">'+th+'<tbody>'+rws+'</tbody></table>')
VALJ=json.dumps({'ours':val['ours'],'bc':val['bc'],'plate':plate})
for k,v in [('/*CSS*/',css),('/*VAL*/',VALJ),('/*CHARTS*/',charts)]:
  assert t.count(k)==1,k;assert '</script' not in v.lower();t=t.replace(k,v)
art=t.replace('/*TOOLURL*/','https://bymankirat.github.io/swapseq/')
open(B+'validation.html','w').write(art);print('artifact',len(art))
gh=t.replace('/*TOOLURL*/','./')
gh=re.sub(r'<link rel="preconnect"[^>]*>(<link rel="preconnect"[^>]*>)?\n?','',gh)
gh=re.sub(r'<link rel="stylesheet" href="https://fonts.googleapis.com[^>]*>','',gh)
ff='''<style>
@font-face{font-family:"Bricolage Grotesque";src:url("fonts/BricolageGrotesque.ttf") format("truetype");font-weight:200 800;font-stretch:75% 100%;font-display:swap}
@font-face{font-family:"Instrument Sans";src:url("fonts/InstrumentSans.ttf") format("truetype");font-weight:400 700;font-stretch:75% 100%;font-display:swap}
@font-face{font-family:"Instrument Sans";src:url("fonts/InstrumentSans-Italic.ttf") format("truetype");font-weight:400 700;font-style:italic;font-display:swap}
@font-face{font-family:"JetBrains Mono";src:url("fonts/JetBrainsMono.ttf") format("truetype");font-weight:100 800;font-display:swap}
</style>'''
gh=gh.replace('<style>',ff+'\n<style>',1)
assert 'googleapis' not in gh and 'gstatic' not in gh
head='<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">\n'
i=gh.index('<style>');j=gh.index('</style>',gh.index('<style>',i+7))+8
gh=head+gh[:j]+'\n</head>\n<body>\n'+gh[j:]+'\n</body>\n</html>\n'
open(B+'validation_gh.html','w').write(gh);print('gh',len(gh))
