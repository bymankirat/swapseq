import json,glob,os,sys
# Usage: python3 atlas/code/val11_sum.py <folder with val11 result .json files> [output .json]
IN=sys.argv[1] if len(sys.argv)>1 else 'val11'
OUT=sys.argv[2] if len(sys.argv)>2 else 'val11_summary.json'
R=[json.load(open(f)) for f in glob.glob(os.path.join(IN,'*.json'))]
def mae(pairs):
    p=[(a,b) for a,b in pairs if a is not None and b is not None];return (sum(abs(a-b) for a,b in p)/len(p)*100, sum(abs(a-b)>.25 for a,b in p), len(p)) if p else None
c_old=mae([(r['clean']['old'],r['clean']['truth']) for r in R]);c_new=mae([(r['clean']['new'],r['clean']['truth']) for r in R])
e_old=mae([(r['withErrors']['old'],r['withErrors']['actual']) for r in R if r['withErrors']]);e_new=mae([(r['withErrors']['new'],r['withErrors']['actual']) for r in R if r['withErrors']])
# tables that are checkable (truth>=0.5) shown as near-zero
miss_old=sum(1 for r in R if r['clean']['truth'] and r['clean']['truth']>=.5 and (r['clean']['old'] or 0)<.2)
miss_new=sum(1 for r in R if r['clean']['truth'] and r['clean']['truth']>=.5 and r['clean']['new'] is not None and r['clean']['new']<.2)
over_old=sum(1 for r in R if r['clean']['truth'] is not None and r['clean']['truth']<.2 and (r['clean']['old'] or 0)>=.5)
over_new=sum(1 for r in R if r['clean']['truth'] is not None and r['clean']['truth']<.2 and (r['clean']['new'] or 0)>=.5)
secs=sorted((r['powerSeconds'],r['n'],r['key']) for r in R)
print('tables',len(R))
print('clean labels  MAE old %.1f (misses>25pt: %d/%d) | new %.1f (%d/%d)'%(c_old[0],c_old[1],c_old[2],c_new[0],c_new[1],c_new[2]))
print('1%% errors     MAE old %.1f (%d/%d) | new %.1f (%d/%d)'%(e_old[0],e_old[1],e_old[2],e_new[0],e_new[1],e_new[2]))
print('checkable tables shown <20%%: old',miss_old,'new',miss_new,'| uncheckable shown >=50%%: old',over_old,'new',over_new)
print('extra seconds: median',secs[len(secs)//2],'max',secs[-1])
json.dump({'clean':{'old':c_old,'new':c_new},'errors':{'old':e_old,'new':e_new},'missOld':miss_old,'missNew':miss_new,'overOld':over_old,'overNew':over_new,'rows':R},open(OUT,'w'),indent=1)
