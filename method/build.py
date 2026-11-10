# Rebuilds the two generated files from method/core.mjs, method/bridge.js and method/app.js:
#   method/swapseq-core.js  (engine + bridge, used by the R package and nothing else to edit by hand)
#   index.html              (the engine and the page script are spliced into the published page)
# Run from the repository root:  python3 method/build.py
import re,sys
core=open('method/core.mjs').read();bridge=open('method/bridge.js').read()
flat=re.sub(r'^export ','',core,flags=re.M)
open('method/swapseq-core.js','w').write('// SPDX-License-Identifier: GPL-3.0-or-later\n// SwapSeq core 1.1.1 (method v3) bundled for R and the command line. Generated from core.mjs; do not edit by hand.\n'+flat+('' if flat.endswith('\n') else '\n')+bridge)
page=open('index.html').read()
a='<script type="text/plain" id="coreSrc">';i=page.index(a)+len(a);j=page.index('</script>',i)
assert '</script' not in flat.lower()
page=page[:i]+flat+page[j:]
app=open('method/app.js').read();assert '</script' not in app.lower()
k=page.index(app[:120]) if app[:120] in page else page.index('// ---------- state ----------')-len('// SPDX-License-Identifier: GPL-3.0-or-later\n// Copyright (C) 2026 Mankirat Singh\n')
e=page.rindex('</script>')
page=page[:k]+app+('' if app.endswith('\n') else '\n')+page[e:]
open('index.html','w').write(page);print('built: swapseq-core.js',len(flat)+len(bridge),'index.html',len(page))
