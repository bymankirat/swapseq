# Build a standalone SwapSeq Atlas page for GitHub Pages from the dashboard page.
import json,re,csv,io
A='/tmp/claude-0/-home-claude/d674b933-5c36-5fd7-9c63-4fe2a591bf95/scratchpad/atlas/'
page=open(A+'index.html').read()
rows=json.load(open(A+'atlas_rows_v2.json'))
page=page.replace('<meta name="dash-one-writer">\n','')
# standalone wording for the data card
old='<p>Every figure on this page comes from one table of study-level results. Open Sources at the top of the page to see it and download it as CSV. Sample IDs are not included; they go to each study\'s authors first.</p>'
new='<p>Every figure on this page comes from one table of study-level results: <a href="data/atlas.csv">atlas.csv</a> · <a href="data/atlas.json">atlas.json</a>. Code: <a href="https://github.com/bymankirat/swapseq/tree/main/atlas">github.com/bymankirat/swapseq/atlas</a>. Sample IDs are not included; they go to each study\'s authors first.</p>'
assert page.count(old)==1;page=page.replace(old,new)
old='<p>Check your own study in the browser at bymankirat.github.io/swapseq, or in R with the swapseqr package.</p>'
new='<p>Check your own study in the browser with <a href="https://bymankirat.github.io/swapseq/">SwapSeq</a>, or in R with the <a href="https://github.com/bymankirat/swapseqr">swapseqr</a> package.</p>'
assert page.count(old)==1;page=page.replace(old,new)
tokens='''
:root{
  --font-anthropic-sans:"IBM Plex Sans",system-ui,-apple-system,"Segoe UI",sans-serif;
  --font-anthropic-serif:"Newsreader",Georgia,"Times New Roman",serif;
  --color-bg:#f6f7f4; --color-panel:#ffffff; --color-fg:#1b1f1c; --color-fg-muted:#5d665f; --color-border-line:#dfe3dc;
  --color-ok:#2f7d4f; --color-warn:#a8650c; --color-bad:#b3362c;
  --cds-chart-grid:#e7eae4; --cds-chart-axis:#9aa39b; --cds-chart-muted:#b9c1ba; --cds-chart-reference:#7d877f;
  --cds-surface-popover:#ffffff; --cds-border:#d5dad2; --cds-shadow-popover:0 6px 24px rgba(20,30,22,.12); --cds-text-secondary:#5d665f;
  --c1:#2f6fd0; --c2:#e07b39; --c3:#2a9d77; --c4:#c2477e; --c5:#8a63c8; --c6:#d1a21f; --c7:#3a9bb5; --c8:#7a7f2a;
  --cds-font-size-caption:.8125rem; --cds-font-size-body:1rem; --cds-font-size-heading:1.125rem; --cds-font-size-title:1.9rem;
  --cds-font-weight-medium:600; --cds-gap-xs:.25rem; --cds-gap-sm:.5rem; --cds-gap-md:1rem; --cds-gap-lg:1.5rem; --cds-gap-xl:2rem;
  --cds-pad-xs:.25rem; --cds-pad-sm:.5rem; --cds-pad-md:.875rem; --cds-pad-lg:1.25rem; --cds-radius:10px; --cds-dur-slow:.6s; --cds-ease-out:cubic-bezier(.2,.8,.2,1);
}
@media (prefers-color-scheme: dark){ :root:not([data-theme="light"]){
  --color-bg:#111412; --color-panel:#181c19; --color-fg:#e8ece7; --color-fg-muted:#9aa49c; --color-border-line:#2b312c;
  --color-ok:#6cc48f; --color-warn:#e3a54f; --color-bad:#ef7b70;
  --cds-chart-grid:#252b26; --cds-chart-axis:#566058; --cds-chart-muted:#4a544c; --cds-chart-reference:#7f8a81;
  --cds-surface-popover:#1d221e; --cds-border:#333a34; --cds-shadow-popover:0 6px 24px rgba(0,0,0,.5); --cds-text-secondary:#9aa49c;
  --c1:#6ea2ff; --c2:#ff9d5c; --c3:#4cc79c; --c4:#ef78ae; --c5:#b394f0; --c6:#f0c64a; --c7:#5cc3dc; --c8:#b6bb55; color-scheme:dark } }
:root[data-theme="dark"]{
  --color-bg:#111412; --color-panel:#181c19; --color-fg:#e8ece7; --color-fg-muted:#9aa49c; --color-border-line:#2b312c;
  --color-ok:#6cc48f; --color-warn:#e3a54f; --color-bad:#ef7b70;
  --cds-chart-grid:#252b26; --cds-chart-axis:#566058; --cds-chart-muted:#4a544c; --cds-chart-reference:#7f8a81;
  --cds-surface-popover:#1d221e; --cds-border:#333a34; --cds-shadow-popover:0 6px 24px rgba(0,0,0,.5); --cds-text-secondary:#9aa49c;
  --c1:#6ea2ff; --c2:#ff9d5c; --c3:#4cc79c; --c4:#ef78ae; --c5:#b394f0; --c6:#f0c64a; --c7:#5cc3dc; --c8:#b6bb55; color-scheme:dark }
*{box-sizing:border-box}
html{scroll-behavior:smooth}
body{margin:0;background:var(--color-bg);color:var(--color-fg);font-family:var(--font-anthropic-sans);font-size:15px;line-height:1.5}
a{color:inherit;text-underline-offset:3px}
#dash-root{max-width:1180px;margin:0 auto;padding-block:20px 60px;padding-left:max(16px,env(safe-area-inset-left,0px));padding-right:max(16px,env(safe-area-inset-right,0px))}
.top{display:flex;justify-content:space-between;align-items:center;gap:12px;flex-wrap:wrap;margin-bottom:16px;font-size:.8125rem;color:var(--color-fg-muted)}
.top b{color:var(--color-fg);font-family:var(--font-anthropic-serif);font-size:1.05rem;font-weight:500}
.top nav{display:flex;gap:14px;flex-wrap:wrap}
.foot{margin-top:40px;padding-top:16px;border-top:1px solid var(--color-border-line);font-size:.8125rem;color:var(--color-fg-muted);display:flex;gap:16px;flex-wrap:wrap;justify-content:space-between}
.dash-select{font:inherit;font-size:.8125rem;padding:.3em 1.6em .3em .7em;border:1px solid var(--color-border-line);border-radius:var(--cds-radius);background-color:var(--color-panel);color:var(--color-fg)}
.sa-search{min-width:0}
'''
shim='''
window.dash=(function(){
  const DATA={atlas:{status:'ok',data:ROWS,columns:Object.keys(ROWS[0]),meta:[]}};
  const calcs={},draws=new Map();
  function data(id){if(DATA[id])return DATA[id];const c=calcs[id];if(c){if(!c.out){try{c.out={status:'ok',data:c.fn(...c.inputs.map(i=>data(i).data)),meta:[]}}catch(e){console.error(e);c.out={status:'error',data:[]}}}return c.out}return {status:'missing',data:[]}}
  function run(name){const d=draws.get(name);if(!d)return;document.querySelectorAll('[data-draw="'+name+'"]').forEach(el=>{if(el.closest('[hidden]'))return;try{if(typeof d.fn==='function')d.fn(el);else{let m=d.made.get(el);if(m===undefined){m=d.fn.mount(el);d.made.set(el,m)}if(d.fn.update)d.fn.update(el,m)}}catch(e){console.error(name,e)}})}
  return {data,calc:(id,o)=>{calcs[id]=o},draw:(name,fn)=>{if(fn===undefined){run(name);return}draws.set(name,{fn,made:new Map()})},
    onData:f=>f(),colors:['var(--c1)','var(--c2)','var(--c3)','var(--c4)','var(--c5)','var(--c6)','var(--c7)','var(--c8)'],
    setLink:id=>{try{history.replaceState(null,'','#'+id)}catch(e){}},params:()=>({}),setParams(){},resetParams(){},refresh(){},loader(){},
    _runAll(){for(const k of draws.keys())run(k)}};
})();
'''
after='''
<script>
  dash._runAll();
  let lastW=innerWidth,rt;addEventListener('resize',()=>{clearTimeout(rt);rt=setTimeout(()=>{if(Math.abs(innerWidth-lastW)>4){lastW=innerWidth;dash._runAll()}},200)});
  matchMedia('(prefers-color-scheme: dark)').addEventListener('change',()=>dash._runAll());
  const go=document.getElementById(location.hash.slice(1));if(go)setTimeout(()=>go.scrollIntoView(),50);
</script>'''
desc='How checkable are public longitudinal microbiome datasets for swapped samples? SwapSeq tested 30 public tables, 16,336 samples, with hidden test swaps.'
html=f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>SwapSeq Atlas</title>
<meta name="description" content="{desc}">
<meta property="og:type" content="website">
<meta property="og:title" content="SwapSeq Atlas: can public microbiome data be checked for swapped samples?">
<meta property="og:description" content="{desc}">
<meta property="og:url" content="https://bymankirat.github.io/swapseq/atlas/">
<meta property="og:image" content="https://bymankirat.github.io/swapseq/atlas/og.png">
<meta property="og:image:width" content="1200"><meta property="og:image:height" content="627">
<meta name="twitter:card" content="summary_large_image">
<meta name="author" content="Mankirat Singh">
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=Newsreader:opsz,wght@6..72,400;6..72,500&display=swap" rel="stylesheet">
<style>{tokens}</style>
<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.8.5/d3.min.js"></script>
</head>
<body>
<div id="dash-root">
<div class="top"><b>SwapSeq Atlas</b><nav><a href="https://bymankirat.github.io/swapseq/">SwapSeq tool</a><a href="https://bymankirat.github.io/swapseq/validation.html">Methods</a><a href="https://github.com/bymankirat/swapseq/tree/main/atlas">Code and data</a></nav></div>
<script>const ROWS={json.dumps(rows,separators=(',',':'))};{shim}</script>
{page}
<div class="foot"><span>Mankirat Singh · <a href="https://orcid.org/0009-0005-3375-2123">ORCID</a> · <a href="https://www.linkedin.com/in/mankirat-singh-694440177">LinkedIn</a></span><span>SwapSeq is free software (GPL-3.0). Public data only; no sample IDs on this page.</span></div>
</div>
{after}
</body>
</html>
'''
open(A+'site/index.html','w').write(html)
import os;os.makedirs(A+'site/data',exist_ok=True)
pub=[{k:v for k,v in r.items()} for r in rows]
json.dump(pub,open(A+'site/data/atlas.json','w'),indent=1)
keys=list(pub[0].keys());buf=io.StringIO();w=csv.DictWriter(buf,fieldnames=keys);w.writeheader();w.writerows(pub);open(A+'site/data/atlas.csv','w').write(buf.getvalue())
print(len(html))
