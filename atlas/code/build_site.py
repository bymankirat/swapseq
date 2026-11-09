# Build a standalone SwapSeq Atlas page for GitHub Pages from the dashboard page.
import json,re,csv,io
A='/tmp/claude-0/-home-claude/d674b933-5c36-5fd7-9c63-4fe2a591bf95/scratchpad/atlas/'
page=open(A+'index.html').read()
rows=json.load(open(A+'atlas_rows_v3.json'))
page=page.replace('<meta name="dash-one-writer">\n','')
page=page.replace('var(--font-anthropic-sans)','var(--sans)').replace('var(--font-anthropic-serif)','var(--display)')
_a="const ink = SA.rgb('var(--color-fg)');"
assert page.count(_a)==1
page=page.replace(_a,"const ink = getComputedStyle(el).color;")
# standalone wording for the data card
old='<p>Every figure on this page comes from one table of study-level results. Open Sources at the top of the page to see it and download it as CSV. Sample IDs are not included; they go to each study\'s authors first.</p>'
new='<p>Every figure on this page comes from one table of study-level results: <a href="data/atlas.csv">atlas.csv</a> · <a href="data/atlas.json">atlas.json</a>. Code: <a href="https://github.com/bymankirat/swapseq/tree/main/atlas">github.com/bymankirat/swapseq/atlas</a>. Sample IDs are not included; they go to each study\'s authors first.</p>'
assert page.count(old)==1;page=page.replace(old,new)
old='<p>Check your own study in the browser at bymankirat.github.io/swapseq, or in R with the swapseqr package.</p>'
new='<p>Check your own study in the browser with <a href="https://bymankirat.github.io/swapseq/">SwapSeq</a>, or in R with the <a href="https://github.com/bymankirat/swapseqr">swapseqr</a> package.</p>'
assert page.count(old)==1;page=page.replace(old,new)
tokens='''
@font-face{font-family:"Bricolage Grotesque";src:url("../fonts/BricolageGrotesque.ttf") format("truetype");font-weight:200 800;font-stretch:75% 100%;font-display:swap}
@font-face{font-family:"Instrument Sans";src:url("../fonts/InstrumentSans.ttf") format("truetype");font-weight:400 700;font-stretch:75% 100%;font-display:swap}
@font-face{font-family:"JetBrains Mono";src:url("../fonts/JetBrainsMono.ttf") format("truetype");font-weight:100 800;font-display:swap}
:root{
  --display:"Bricolage Grotesque","Avenir Next","Segoe UI",system-ui,sans-serif;
  --sans:"Instrument Sans","Segoe UI",system-ui,-apple-system,sans-serif;
  --mono:"JetBrains Mono",ui-monospace,SFMono-Regular,Menlo,monospace;
  --color-bg:#F3F5F8; --color-panel:#FFFFFF; --color-fg:#172036; --color-fg-muted:#56617A; --color-border-line:#DCE1E9;
  --color-ok:#1E8458; --color-warn:#A86200; --color-bad:#D2304A; --accent:#2B47C9; --accent-ink:#FFFFFF;
  --cds-chart-grid:#E4E8EF; --cds-chart-axis:#9AA3B5; --cds-chart-muted:#B9C1CF; --cds-chart-reference:#7C869B;
  --cds-surface-popover:#FFFFFF; --cds-border:#D3D9E3; --cds-shadow-popover:0 8px 28px rgba(20,30,60,.14); --cds-text-secondary:#56617A;
  --c1:#2B47C9; --c2:#E5762F; --c3:#1E9E74; --c4:#D2304A; --c5:#7B5BD6; --c6:#D9A21B; --c7:#1F97B8; --c8:#7C8A25;
  --cds-font-size-caption:.8125rem; --cds-font-size-body:1rem; --cds-font-size-heading:1.125rem; --cds-font-size-title:2rem;
  --cds-font-weight-medium:600; --cds-gap-xs:.25rem; --cds-gap-sm:.5rem; --cds-gap-md:1rem; --cds-gap-lg:1.5rem; --cds-gap-xl:2rem;
  --cds-pad-xs:.25rem; --cds-pad-sm:.5rem; --cds-pad-md:.875rem; --cds-pad-lg:1.25rem; --cds-radius:12px; --cds-dur-slow:.6s; --cds-ease-out:cubic-bezier(.2,.8,.2,1);
}
@media (prefers-color-scheme: dark){ :root:not([data-theme="light"]){
  --color-bg:#0E1424; --color-panel:#151D31; --color-fg:#E9EDF6; --color-fg-muted:#A0AAC0; --color-border-line:#263049;
  --color-ok:#56D195; --color-warn:#F3B24C; --color-bad:#FF6F80; --accent:#93A4FF; --accent-ink:#0E1424;
  --cds-chart-grid:#222C44; --cds-chart-axis:#4A5677; --cds-chart-muted:#3A4664; --cds-chart-reference:#7C869B;
  --cds-surface-popover:#1A2338; --cds-border:#2E3955; --cds-shadow-popover:0 8px 28px rgba(0,0,0,.5); --cds-text-secondary:#A0AAC0;
  --c1:#93A4FF; --c2:#FF9D5C; --c3:#4CD3A2; --c4:#FF6F80; --c5:#B79CFF; --c6:#F3C64A; --c7:#5CC8E6; --c8:#C2CF55; color-scheme:dark } }
:root[data-theme="dark"]{
  --color-bg:#0E1424; --color-panel:#151D31; --color-fg:#E9EDF6; --color-fg-muted:#A0AAC0; --color-border-line:#263049;
  --color-ok:#56D195; --color-warn:#F3B24C; --color-bad:#FF6F80; --accent:#93A4FF; --accent-ink:#0E1424;
  --cds-chart-grid:#222C44; --cds-chart-axis:#4A5677; --cds-chart-muted:#3A4664; --cds-chart-reference:#7C869B;
  --cds-surface-popover:#1A2338; --cds-border:#2E3955; --cds-shadow-popover:0 8px 28px rgba(0,0,0,.5); --cds-text-secondary:#A0AAC0;
  --c1:#93A4FF; --c2:#FF9D5C; --c3:#4CD3A2; --c4:#FF6F80; --c5:#B79CFF; --c6:#F3C64A; --c7:#5CC8E6; --c8:#C2CF55; color-scheme:dark }
*{box-sizing:border-box}
html{scroll-behavior:smooth}
body{margin:0;background:var(--color-bg);color:var(--color-fg);font:16px/1.55 var(--sans);-webkit-font-smoothing:antialiased}
a{color:var(--accent);text-underline-offset:3px}
.ss-top{position:sticky;top:0;z-index:20;background:color-mix(in srgb,var(--color-bg) 88%,transparent);backdrop-filter:blur(10px);-webkit-backdrop-filter:blur(10px);border-bottom:1px solid var(--color-border-line)}
.ss-top .in{max-width:1180px;margin:0 auto;display:flex;align-items:center;gap:24px;padding-block:12px;padding-left:max(16px,env(safe-area-inset-left,0px));padding-right:max(16px,env(safe-area-inset-right,0px))}
.ss-brand{display:flex;align-items:center;gap:10px;font:700 20px/1 var(--display);letter-spacing:-.02em;color:var(--color-fg);text-decoration:none}
.ss-brand small{font:500 13px var(--sans);color:var(--color-fg-muted);letter-spacing:0}
.ss-top nav{display:flex;gap:18px;margin-left:auto;flex-wrap:wrap}
.ss-top nav a{color:var(--color-fg-muted);text-decoration:none;font-size:15px;white-space:nowrap}
.ss-top nav a:hover,.ss-top nav a[aria-current="page"]{color:var(--color-fg)}
@media (max-width:720px){.ss-top nav a:not(.keep){display:none}}
#dash-root{max-width:1180px;margin:0 auto;padding-block:28px 64px;padding-left:max(16px,env(safe-area-inset-left,0px));padding-right:max(16px,env(safe-area-inset-right,0px))}
/* SwapSeq look for the atlas page */
.sa-title,.sa-h2,.sa-profile h3{font-weight:700 !important;letter-spacing:-.025em}
.sa-title{font-size:clamp(2.2rem,5vw,3.6rem) !important;line-height:1.02 !important}
.sa-h2{font-size:clamp(1.6rem,3.2vw,2.2rem) !important}
.sa-tile-num{font-weight:700 !important;letter-spacing:-.02em}
.sa-hero{--color-panel:#0B1122;--color-fg:#F4F6FB;--color-fg-muted:#AAB3CB;--color-border-line:#26314F;color:var(--color-fg);
  background-color:#0B1122 !important;background-image:radial-gradient(rgba(147,164,255,.16) 1px,transparent 1.3px),radial-gradient(120% 90% at 85% 10%,rgba(43,71,201,.35),transparent 60%) !important;
  background-size:20px 20px,100% 100% !important;border-color:#1D2745 !important;min-height:28rem !important}
.sa-hero .sa-kicker{color:#93A4FF}
.sa-hero .sa-lede b{color:#FFFFFF}
.sa-nav{top:57px !important}
.sa-chip[aria-pressed="true"]{background:var(--accent) !important;border-color:var(--accent) !important;color:var(--accent-ink) !important}
.sa-btn{border-color:var(--accent) !important;color:var(--accent) !important}
.foot{margin-top:48px;padding-top:18px;border-top:1px solid var(--color-border-line);font-size:.8125rem;color:var(--color-fg-muted);display:grid;gap:10px;line-height:1.55}
.foot b{color:var(--color-fg)}
.dash-select{font:inherit;font-size:.8125rem;padding:.35em 1.6em .35em .7em;border:1px solid var(--color-border-line);border-radius:10px;background-color:var(--color-panel);color:var(--color-fg)}
.sa-search{min-width:0}
.sa-label svg text{font-size:.9rem !important;font-weight:700}
.sa-section{scroll-margin-top:118px}
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
desc='How checkable are public longitudinal microbiome datasets for swapped samples? SwapSeq tested 36 public tables, 18,317 samples, with hidden test swaps.'
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
<style>{tokens}</style>
<link rel="icon" href="../assets/favicon.svg" type="image/svg+xml">
<script src="vendor/d3.min.js"></script>
</head>
<body>
<header class="ss-top"><div class="in"><a class="ss-brand" href="../"><svg width="28" height="28" viewBox="0 0 28 28" aria-hidden="true"><rect x="1.5" y="1.5" width="25" height="25" rx="6" fill="none" stroke="currentColor" stroke-width="1.6" opacity=".35"/><circle cx="9" cy="14" r="4.6" fill="var(--accent)"/><circle cx="19" cy="14" r="4.6" fill="var(--color-bad)"/><path d="M8 7.2c3-3 9-3 12 0M20 20.8c-3 3-9 3-12 0" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round"/></svg><span>SwapSeq <small>Atlas</small></span></a>
<nav aria-label="Site"><a class="keep" href="../">Tool</a><a href="../validation.html">Methods</a><a class="keep" href="./" aria-current="page">Atlas</a><a href="https://github.com/bymankirat/swapseqr">R package</a><a href="https://github.com/bymankirat/swapseq/tree/main/atlas">Code and data</a></nav></div></header>
<div id="dash-root">
<script>const ROWS={json.dumps(rows,separators=(',',':'))};{shim}</script>
{page}
<div class="foot">
<p style="margin:0"><b>Data.</b> Built only from datasets that their authors made public. Each dataset stays under its own terms; this page shows study-level statistics derived from them, with full credit in the Method section. No sample IDs, no personal data. A flag is a candidate for a records check, not a finding of error; sample-level details go to each study's authors first.</p>
<p style="margin:0"><b>Privacy.</b> No cookies, no analytics, no tracking. Fonts and scripts are served from this site; nothing is loaded from third parties. The website host may keep standard access logs for security.</p>
<p style="margin:0"><b>Use.</b> Research quality control only. Not a medical device; not for diagnosis or clinical decisions. Provided without warranty.</p>
<p style="margin:0"><b>Licences.</b> Code: GPL-3.0-or-later. Fonts: SIL Open Font License (Bricolage Grotesque, Instrument Sans, JetBrains Mono). d3: ISC licence.</p>
<p style="margin:0">Made by Mankirat Singh · <a href="https://orcid.org/0009-0005-3375-2123">ORCID</a> · <a href="https://www.linkedin.com/in/mankirat-singh-694440177">LinkedIn</a> · <a href="https://github.com/bymankirat/swapseq/issues/new/choose">Feedback</a></p>
</div>
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
