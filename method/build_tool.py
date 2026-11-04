import json,re
B='/tmp/claude-0/-home-claude/d674b933-5c36-5fd7-9c63-4fe2a591bf95/scratchpad/swap/'
t=open(B+'tool.html').read()
core=re.sub(r'^export ','',open(B+'core.mjs').read(),flags=re.M)
val={'ours':json.load(open(B+'bench_results_prior0.01_v3.json'))+json.load(open(B+'bench_results_prior0.01_v3new.json'))+json.load(open(B+'bench_results_prior0.01_imp.json')),'bc':json.load(open(B+'bench_results_prior0.01_braycurtis_v3.json'))+json.load(open(B+'bench_results_prior0.01_braycurtis_v3new.json'))+json.load(open(B+'bench_results_prior0.01_braycurtis_imp.json'))}
D=json.load(open(B+'demo.json'))
# simulated plate layout for the demo: samples plated in week order, row by row; the planted swap pair sits side by side
rows=[l.split('\t') for l in D['meta'].strip().split('\n')];hdr,body=rows[0],rows[1:]
a,b=D['answer']['swap'];body.sort(key=lambda r:(float(r[2]),r[1]))
ia=[r[0] for r in body].index(a);rb=[r for r in body if r[0]==b][0];body.remove(rb);ia=[r[0] for r in body].index(a);body.insert(ia+1,rb)
out=['\t'.join(hdr+['plate','well'])]
for k,r in enumerate(body):out.append('\t'.join(r+['plate%d'%(k//96+1),'%s%02d'%('ABCDEFGH'[(k%96)//12],k%12+1)]))
wa=out[ia+1].split('\t')[-1];wb=out[ia+2].split('\t')[-1];assert wa[0]==wb[0],(wa,wb)
D['meta']='\n'.join(out)+'\n';demo=json.dumps(D)

# ---------- hero: the demo plate, with row E opened up as a rack of tubes ----------
# cap colour = the name written on the tube; liquid colour = whose microbes are inside
def plate_svg():
  rs=[l.split('\t') for l in D['meta'].strip().split('\n')][1:]
  lab={r[0]:r[1] for r in rs}
  plate=[r for r in rs if r[4]=='plate1']
  row=sorted([r for r in plate if r[5].startswith('E')],key=lambda r:r[5])[:8]
  a,b=D['answer']['swap'];c=D['answer']['relabel']
  true={a:lab[b],b:lab[a],c:'P03'}
  pal=['hsl(210 62% 54%)','hsl(30 82% 56%)','hsl(266 46% 62%)','hsl(150 44% 44%)','hsl(186 52% 44%)','hsl(332 52% 62%)','hsl(46 74% 50%)','hsl(6 56% 56%)','hsl(96 38% 48%)']
  people=[]
  for r in row:
    for p_ in (r[1],true.get(r[0],r[1])):
      if p_ not in people:people.append(p_)
  col={p_:pal[i%len(pal)] for i,p_ in enumerate(people)}
  rest=sorted(set(r[1] for r in rs)-set(col))
  for i,p_ in enumerate(rest):col[p_]='hsl(%d 20%% 64%%)'%((i*137.5+80)%360)
  o=['<svg viewBox="0 0 640 712" role="img" aria-label="The demo study\'s 96-well plate. Row E is opened up as eight tubes: E01 and E02 hold each other\'s samples, and E03 holds a sample from a different person than its label says.">','<defs>']
  X=lambda i:46+i*78;OY=338
  for i in range(len(row)):o.append('<clipPath id="tb%d"><path d="M%d %d V%d a20 20 0 0 0 40 0 V%d Z"/></clipPath>'%(i,X(i)-20,OY+96,OY+262,OY+96))
  o.append('</defs>')
  # plate
  px0,pd,py0=133,34,62
  # zoom connectors from row E wells down to the tubes
  o.append('<g class="zoom" fill="none" stroke="var(--accent)" stroke-width="1.2" stroke-opacity=".45">')
  for i in range(len(row)):o.append('<path d="M%d %d C%d %d %d %d %d %d"/>'%(px0+i*pd,py0+4*pd+13,px0+i*pd,py0+4*pd+80,X(i),OY+10,X(i),OY+60))
  o.append('</g>')
  o.append('<path d="M126 6 H560 Q574 6 574 20 V300 Q574 314 560 314 H80 Q66 314 66 300 V66 Z" fill="var(--plate)" stroke="var(--line2)" stroke-width="1.5"/>')
  for k in range(12):o.append('<text x="%d" y="38" text-anchor="middle" font-size="11" fill="var(--faint)" font-family="var(--sans)">%d</text>'%(px0+k*pd,k+1))
  for r_ in range(8):o.append('<text x="96" y="%d" text-anchor="middle" font-size="11" fill="var(--faint)" font-family="var(--sans)">%s</text>'%(py0+r_*pd+4,'ABCDEFGH'[r_]))
  o.append('<rect class="scanbox" x="%d" y="%d" width="%d" height="%d" rx="15" fill="var(--accent-soft)" stroke="var(--accent)" stroke-width="2"/>'%(px0-22,py0-16,11*pd+44,32))
  for r in plate:
    w=r[5];ri='ABCDEFGH'.index(w[0]);ci=int(w[1:])-1;cx,cy=px0+ci*pd,py0+ri*pd
    inrow=(ri==4 and ci<8);bad=r[0] in true;kind='swap' if r[0] in (a,b) else 'mis'
    o.append('<circle class="pw%s" data-r="%d" data-c="%d" cx="%d" cy="%d" r="12" fill="%s" stroke="var(--rim)" stroke-width="1"><title>%s, labelled %s</title></circle>'%(' inrow' if inrow else '',ri,ci,cx,cy,col[r[1]],w,r[1]))
    if bad:o.append('<circle class="pring" cx="%d" cy="%d" r="15.5" fill="none" stroke="var(--%s)" stroke-width="2.6"/>'%(cx,cy,'swap' if kind=='swap' else 'mis'))
  o.append('<line x1="8" x2="632" y1="%d" y2="%d" stroke="var(--line2)" stroke-width="1.5"/>'%(OY+301,OY+301))
  for i,r in enumerate(row):
    sid,w,l=r[0],r[5],r[1];t=true.get(sid,l);cx=X(i);bad=sid in true
    kind='swap' if sid in (a,b) else ('mis' if bad else 'ok')
    o.append('<g class="tube">')
    if bad:o.append('<rect class="flagbox" x="%d" y="%d" width="58" height="236" rx="16" fill="none" stroke="var(--%s)" stroke-width="2.5"/>'%(cx-29,OY+62,'swap' if kind=='swap' else 'mis'))
    o.append('<path d="M%d %d V%d a20 20 0 0 0 40 0 V%d Z" fill="var(--plate)" stroke="var(--line2)" stroke-width="1.6"/>'%(cx-20,OY+96,OY+262,OY+96))
    o.append('<g clip-path="url(#tb%d)"><rect class="liq" x="%d" y="%d" width="40" height="110" fill="%s"><title>Microbes inside: %s</title></rect></g>'%(i,cx-20,OY+176,col[t],t))
    o.append('<rect x="%d" y="%d" width="48" height="28" rx="6" fill="%s"><title>Name on the tube: %s</title></rect><rect x="%d" y="%d" width="52" height="8" rx="3" fill="%s"/>'%(cx-24,OY+70,col[l],l,cx-26,OY+66,col[l]))
    o.append('<rect x="%d" y="%d" width="36" height="24" rx="4" fill="var(--surface)" stroke="var(--line2)"/><text x="%d" y="%d" text-anchor="middle" font-size="13" font-weight="600" fill="var(--ink)" font-family="var(--mono)">%s</text>'%(cx-18,OY+112,cx,OY+129,l))
    o.append('<text x="%d" y="%d" text-anchor="middle" font-size="14" fill="var(--muted)" font-family="var(--sans)">%s</text>'%(cx,OY+322,w))
    if bad:o.append('<text class="flagtxt" x="%d" y="%d" text-anchor="middle" font-size="13.5" font-weight="700" fill="var(--%s)" font-family="var(--sans)">holds %s</text>'%(cx,OY+348,'swap' if kind=='swap' else 'mis',t))
    o.append('</g>')
  x0,x1=X(0),X(1);m=(x0+x1)/2;Y=OY
  o.append('<g class="swaparc"><path d="M%d %d C%d %d %d %d %d %d" fill="none" stroke="var(--swap)" stroke-width="2.6" stroke-linecap="round"/><path d="M%d %d l-7 -6 M%d %d l2 -9 M%d %d l7 -6 M%d %d l-2 -9" stroke="var(--swap)" stroke-width="2.6" stroke-linecap="round"/><rect x="%d" y="%d" width="84" height="24" rx="12" fill="var(--swap)"/><text x="%d" y="%d" text-anchor="middle" font-size="13.5" font-weight="700" fill="var(--surface)" font-family="var(--sans)">swapped</text></g>'%(x0+6,Y+58,x0+10,Y+22,x1-10,Y+22,x1-6,Y+58,x0+6,Y+58,x0+6,Y+58,x1-6,Y+58,x1-6,Y+58,m-42,Y+2,m,Y+19))
  o.append('</svg>');return '\n'.join(o)

# ---------- validation chart ----------
def val_svg():
  name=lambda d:(d.replace('SMP ','SMP ').replace(' (2 per person)',', 2 per person').replace('IMP adults (16S)','IMP gut, monthly').replace('Diet study (strain-level)','Diet study, strain level').replace('Diet study (species)','Diet study, species').replace('MDSINE2 mice (healthy donor)','Mice, healthy donor').replace('MDSINE2 mice (uc donor)','Mice, UC donor'))
  def grp(d):
    if '2 per person' in d:return 'Before/after, 2 samples per person' if 'Diet' not in d else 'Shotgun, adults sampled daily'
    if d.startswith('SMP') or d.startswith('IMP'):return '16S, adults'
    if d.startswith('Diet'):return 'Shotgun, adults sampled daily'
    return 'Little individual signal (controls)'
  order=['SMP gut','SMP forehead','SMP palm','IMP adults (16S)','SMP tongue','SMP gut (2 per person)','SMP forehead (2 per person)','Diet study (species)','Diet study (strain-level)','Diet study (strain-level) (2 per person)','ECAM infants','MDSINE2 mice (healthy donor)','MDSINE2 mice (uc donor)']
  by={o['dataset']:(o,bc) for o,bc in zip(val['ours'],val['bc'])}
  rows=[(d,)+by[d] for d in order if d in by]
  s16=[o['swappedSamplesCaught'] for d,o,bc in rows if grp(d) in ('16S, adults','Before/after, 2 samples per person')]
  b16=[bc['swappedSamplesCaught'] for d,o,bc in rows if grp(d) in ('16S, adults','Before/after, 2 samples per person')]
  x0,x1,y0,rh,gh=270,736,112,26,34
  X=lambda v:x0+v*(x1-x0)
  ys=[];g=None;y=y0
  for d,o,bc in rows:
    if grp(d)!=g:
      if g is not None:y+=gh
      g=grp(d);ys.append(('h',g,y-22))
    ys.append(('r',(d,o,bc),y));y+=rh
  H=y+40
  O=['<svg viewBox="0 0 760 %d" role="img" aria-label="Share of planted swaps flagged by SwapSeq and by Bray-Curtis in each dataset" font-family="var(--sans)" font-size="13">'%H]
  O.append('<text x="0" y="24" font-family="var(--display)" font-size="19" font-weight="700" fill="var(--ink)">On 16S data SwapSeq flags %d–%d%% of swapped samples, Bray-Curtis %d–%d%%</text>'%(round(min(s16)*100),round(max(s16)*100),round(min(b16)*100),round(max(b16)*100)))
  O.append('<g><circle cx="7" cy="58" r="5.5" fill="var(--obs)"/><text x="19" y="62" fill="var(--muted)">Bray-Curtis, same pipeline</text><circle cx="222" cy="58" r="6.5" fill="var(--accent)"/><text x="235" y="62" fill="var(--muted)">SwapSeq strong flag</text><circle cx="398" cy="58" r="6.5" fill="none" stroke="var(--accent)" stroke-width="2.2"/><text x="411" y="62" fill="var(--muted)">SwapSeq flag or check list</text></g>')
  for t in [0,.25,.5,.75,1]:O.append('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="var(--grid)"/><text x="%.1f" y="%d" text-anchor="middle" fill="var(--faint)">%d%%</text>'%(X(t),X(t),y0-26,y-6,X(t),y+14,round(t*100)))
  for kind,val_,yy in ys:
    if kind=='h':O.append('<text x="0" y="%d" font-weight="700" fill="var(--ink)">%s</text>'%(yy,val_));continue
    d,o,bc=val_;sw,ck,b=o['swappedSamplesCaught'],o['withReview']['swaps'],bc['swappedSamplesCaught'];n=name(d)
    O.append('<g><text x="14" y="%d" fill="var(--ink)">%s</text><line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="var(--line2)" stroke-width="2.5" stroke-linecap="round"/>'%(yy+4,n,X(min(b,sw)),X(max(ck,b)),yy,yy))
    O.append('<circle cx="%.1f" cy="%d" r="5.5" fill="var(--obs)"><title>%s: Bray-Curtis flags %d%%</title></circle>'%(X(b),yy,n,round(b*100)))
    O.append('<circle cx="%.1f" cy="%d" r="6.5" fill="var(--surface)" stroke="var(--accent)" stroke-width="2.2"><title>%s: SwapSeq flag or check list %d%%</title></circle>'%(X(ck),yy,n,round(ck*100)))
    O.append('<circle cx="%.1f" cy="%d" r="6.5" fill="var(--accent)"><title>%s: SwapSeq strong flag %d%%</title></circle></g>'%(X(sw),yy,n,round(sw*100)))
  O.append('</svg>');return '\n'.join(O)

plate={k:json.load(open(B+f'plate_results_burst_{f}.json'))['res'][n] for k,f,n in [('neighbour','L5','neighbour swap'),('L3','L3','shift run'),('L5','L5','shift run'),('L8','L8','shift run')]}
for k in ['L3']:plate[k]['shiftNamed']=plate[k].get('shiftNamed',0)
plate['neighbour']['shiftNamed']=None
VALJ=json.dumps({'ours':val['ours'],'bc':val['bc'],'plate':plate})
for k,v in [('/*PLATE*/',plate_svg()),('/*VALCHART*/',val_svg()),('/*VAL*/',VALJ),('/*JSZIP*/',open(B+'lib/jszip/dist/jszip.min.js').read()),('/*HDF5*/',open(B+'lib/hdf5.min.js').read()),('/*CORE*/',core),('/*DEMO*/',demo),('/*APP*/',open(B+'app.js').read())]:
    assert t.count(k)==1,k;assert '</script' not in v.lower(),k;t=t.replace(k,v)
pc=lambda v:'–' if v is None else '%d%%'%round(v*100)
rws=''
for o,bc in zip(val['ours'],val['bc']):
  w=o['withReview'];rws+=f"<tr><td>{o['dataset']}</td><td class='num'>{o['samples']}</td><td class='num'>{o['people']}</td><td class='num best'>{pc(o['swappedSamplesCaught'])}</td><td class='num'>{pc(w['swaps'])}</td><td class='num'>{pc(bc['swappedSamplesCaught'])}</td><td class='num'>{pc(o['oneWayCaught'])}</td><td class='num'>{pc(w['oneWay'])}</td><td class='num'>{o['falseFlagsPer1000']:.1f}</td><td class='num'>{w['falsePer1000']:.1f}</td></tr>"
th="<thead><tr><th>Dataset</th><th class='num'>Samples</th><th class='num'>People</th><th class='num'>Swaps flagged</th><th class='num'>+ check list</th><th class='num'>Bray-Curtis</th><th class='num'>Mislabels flagged</th><th class='num'>+ check list</th><th class='num'>False flags / 1000</th><th class='num'>False checks / 1000</th></tr></thead>"
old='<table class="val" id="valTable"></table>';assert t.count(old)==1;t=t.replace(old,'<table class="val" id="valTable">'+th+'<tbody>'+rws+'</tbody></table>')
open(B+'swapseq.html','w').write(t);print(len(t))
# GitHub Pages build: same page, fonts served from its own folder instead of Google (no visitor IP goes to Google)
import re as _re
gh=_re.sub(r'<link rel="preconnect"[^>]*>(<link rel="preconnect"[^>]*>)?\n?','',t)
gh=_re.sub(r'<link rel="stylesheet" href="https://fonts.googleapis.com[^>]*>','',gh)
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
open(B+'swapseq_gh.html','w').write(gh);print('gh',len(gh))
