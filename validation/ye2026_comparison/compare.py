# Ye et al. (2026) Find-abnormality vs SwapSeq on the same 15 public metagenomic datasets.
# 1) published labels: which samples each method flags; 2) hidden test swaps: which method finds them.
import json, os, random, subprocess, sys, csv
import os; S=os.environ.get('YE_WORK','.')
D=S+'/ye_data/Find-mislabels_Profiles_and_metadata'; YE=S+'/ye_repo/find_abnormality.py'; OUT=S+'/yecmp/out'
os.makedirs(OUT,exist_ok=True)
REPS=int(os.environ.get('REPS','5'))
def ye(prof,meta,tag):
    wd=f'{OUT}/ye_{tag}'; os.makedirs(wd,exist_ok=True)
    subprocess.run(['python3','-I',YE,prof,meta,'-s',tag],cwd=wd,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=1800)
    mis,dup=set(),set()
    p=f'{wd}/{tag}_true_labels.txt'
    if os.path.exists(p):
        for l in open(p).read().splitlines()[1:]:
            if l.strip(): mis.add(l.split('\t')[0].split(',')[0])
    p=f'{wd}/{tag}_cheating_sns.txt'
    if os.path.exists(p):
        for l in open(p):
            if '==' in l: dup.update(x.strip() for x in l.split('=='))
    return mis,dup
def ss(prof,meta,tag,power):
    o=f'{OUT}/ss_{tag}.json'
    subprocess.run(['node','--max-old-space-size=6000',S+'/yecmp/swapseq_run.mjs',prof,meta,o,'1' if power else '0'],check=True,timeout=7200)
    return json.load(open(o))
def readmeta(p):
    rows=list(csv.reader(open(p),delimiter='\t')); return rows[0],rows[1:]
res=[]
names=sorted(f[:-9] for f in os.listdir(D) if f.endswith('_meta.txt'))
if len(sys.argv)>1: names=[n for n in names if n in sys.argv[1:]]
for name in names:
    prof,meta=f'{D}/{name}_profiles.txt',f'{D}/{name}_meta.txt'
    mis,dup=ye(prof,meta,name+'_orig'); s=ss(prof,meta,name+'_orig',True)
    sflag={r['id'] for r in s['rows'] if r['level']=='flag'}; scheck={r['id'] for r in s['rows'] if r['level']=='check'}
    yeall=mis|dup
    head,rows=readmeta(meta); pat={r[0]:r[1] for r in rows}
    profids=set(open(prof).readline().rstrip('\n').split('\t')[1:])
    cnt={}
    for r in rows:
        if r[0] in profids: cnt[r[1]]=cnt.get(r[1],0)+1
    pool=[r[0] for r in rows if r[0] in profids and cnt[r[1]]>=3 and r[0] not in yeall and r[0] not in sflag and r[0] not in scheck]
    n=s['n']; npairs=int(os.environ['NP']) if os.environ.get('NP') else (1 if n<300 else max(3,round(n*0.003)))
    caught={'ye':0,'ss':0,'ss_check':0,'tot':0}; fp={'ye':0,'ss':0,'ss_check':0,'tot':0}
    for rep in range(REPS):
        R=random.Random(1000+rep); lab=dict(pat); used=set(); planted=[]
        tries=0
        while len(planted)<2*npairs and tries<10000 and len(pool)>=4:
            tries+=1; a,b=R.sample(pool,2)
            if a in used or b in used or lab[a]==lab[b]: continue
            lab[a],lab[b]=lab[b],lab[a]; used|={a,b}; planted+=[a,b]
        if not planted: break
        mp=f'{OUT}/{name}_rep{rep}_meta.txt'
        with open(mp,'w',newline='') as f:
            w=csv.writer(f,delimiter='\t'); w.writerow(head)
            for r in rows: w.writerow([r[0],lab[r[0]]]+r[2:])
        m2,d2=ye(prof,mp,f'{name}_rep{rep}'); s2=ss(prof,mp,f'{name}_rep{rep}',False)
        y2=m2|d2; f2={r['id'] for r in s2['rows'] if r['level']=='flag'}; c2={r['id'] for r in s2['rows'] if r['level']!='ok'}
        for x in planted:
            caught['tot']+=1; caught['ye']+=x in y2; caught['ss']+=x in f2; caught['ss_check']+=x in c2
        others=[x for x in pool if x not in used]
        for x in others:
            fp['tot']+=1; fp['ye']+=x in y2; fp['ss']+=x in f2; fp['ss_check']+=x in c2
    pw=s['power'] or {}
    r=dict(dataset=name,samples=n,people=s['subjects'],
        ye_flagged=len(yeall),ye_mislabel=len(mis),ye_duplicate=len(dup),
        ss_flags=len(sflag),ss_checks=len(scheck),both=len(yeall&sflag),both_incl_check=len(yeall&(sflag|scheck)),
        ss_flag_ids=sorted(sflag),ye_ids=sorted(yeall),
        ss_power_swap=pw.get('swap'),ss_power_check=pw.get('swapWithCheck'),
        hidden_tested=caught['tot'],ye_caught=caught['ye'],ss_caught=caught['ss'],ss_caught_check=caught['ss_check'],
        others_tested=fp['tot'],ye_false=fp['ye'],ss_false=fp['ss'],ss_false_check=fp['ss_check'])
    res.append(r)
    pc=lambda a,b: f'{100*a/b:.0f}%' if b else 'n/a'
    print(f"{name:20s} n={n:5d} | published labels: Ye {len(yeall):3d} (mis {len(mis)}, dup {len(dup)})  SwapSeq {len(sflag):3d} flags +{len(scheck)} checks  overlap {len(yeall&sflag)} | hidden swaps: Ye {pc(caught['ye'],caught['tot'])}  SwapSeq {pc(caught['ss'],caught['tot'])} (+check {pc(caught['ss_check'],caught['tot'])}) | false: Ye {pc(fp['ye'],fp['tot'])} SwapSeq {pc(fp['ss'],fp['tot'])}",flush=True)
    json.dump(res,open(S+'/yecmp/results.json','w'),indent=1)
