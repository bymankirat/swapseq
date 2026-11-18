#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Rerun SwapSeq's main benchmark on your own computer and compare it with the published numbers
# (the README "Validation" table and validation/bench_results_*.json).
#
# Needs: git, unzip, Node.js 18 or newer, and internet for the first run (public data from GitHub,
# only the folders that are needed). Takes about 5 minutes on a laptop.
#
# Usage, from the swapseq folder:
#   bash validation/reproduce.sh              # data go to ~/swapseq-data
#   bash validation/reproduce.sh /some/folder # or choose the data folder
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="${1:-$HOME/swapseq-data}"
OUT="$ROOT/validation/rerun"
mkdir -p "$DATA/data" "$DATA/x" "$OUT"

for tool in git unzip node; do
  command -v "$tool" >/dev/null || { echo "Please install $tool first."; exit 1; }
done
node -e 'process.exit(+process.versions.node.split(".")[0] >= 18 ? 0 : 1)' || { echo "Node.js 18 or newer is needed."; exit 1; }

# get NAME URL COMMIT PATH... : download only the listed paths of a public repository, at a fixed commit
get() {
  local name=$1 url=$2 commit=$3; shift 3
  local d="$DATA/data/$name"
  if [ ! -d "$d/.git" ]; then
    echo "Downloading $name ..."
    git clone -q --filter=blob:none --no-checkout "$url" "$d"
  fi
  git -C "$d" sparse-checkout set --no-cone "$@"
  git -C "$d" checkout -q "$commit"
}

get smp https://github.com/caporaso-lab/student-microbiome-project a60b3c52ac766bb787602b817fac3b3a4dda4208 \
  '/otu-tables/otu_table_even10000_*.ts_only.biom.gz' '/analysis-results/issue-13/map.txt'
get q2-longitudinal https://github.com/qiime2/q2-longitudinal 8c6845e5a71a6c57bce32eecdd1289691d676cbc \
  '/q2_longitudinal/tests/data/ecam_map_maturity.txt' '/q2_longitudinal/tests/data/ecam-table-maturity.qza'
get dietstudy_analyses https://github.com/knights-lab/dietstudy_analyses fc713d5b9d127be61b7d7e6ced55e0506b3a441f \
  '/data/maps/SampleID_map.txt' '/data/microbiome/BURST.tax.txt' '/data/microbiome/processed_sample/taxonomy_counts_s.txt'
get IMP_analyses https://github.com/knights-lab/IMP_analyses 7f5c37092e7bd580e6e6dd8e88db482962090053 \
  '/data/map.txt' '/data/denovo/final_otu.txt'
get MDSINE2_Paper https://github.com/gerberlab/MDSINE2_Paper 1b20972a7dacd9dca9d801d3011556647bbc64f1 \
  '/datasets/gibson/healthy/raw_tables/' '/datasets/gibson/uc/raw_tables/'

if [ ! -d "$DATA/x/ecam-table-maturity" ]; then
  unzip -q "$DATA/data/q2-longitudinal/q2_longitudinal/tests/data/ecam-table-maturity.qza" -d "$DATA/x/ecam-table-maturity"
fi

cd "$OUT"
run() { # run TAG REPS [bc] SETS : one benchmark run, written to validation/rerun/
  local tag=$1 reps=$2 bc=$3 sets=$4
  echo "Running: $sets ($reps repetitions${bc:+, Bray-Curtis})"
  SWAPSEQ_DATA="$DATA" SETS="$sets" TAG="$tag" node "$ROOT/method/bench.mjs" "$reps" 0.01 $bc > /dev/null
}
run _v3      30 ""  gut,tongue,palm,forehead,ecam,gut2,forehead2
run _imp     60 ""  imp
run _v3new   30 ""  diet,dietsp,mdh,mdu,diet2
run _v3      30 bc  gut,tongue,palm,forehead,ecam,gut2,forehead2
run _imp     60 bc  imp

echo
echo "Comparing with the published results:"
node - "$ROOT/validation" "$OUT" <<'JS'
const fs = require('fs'), path = require('path');
const [pub, mine] = process.argv.slice(2);
let same = 0, diff = 0;
for (const f of fs.readdirSync(mine).filter(f => f.endsWith('.json')).sort()) {
  const a = JSON.parse(fs.readFileSync(path.join(pub, f))), b = JSON.parse(fs.readFileSync(path.join(mine, f)));
  for (const r of b) {
    const p = a.find(x => x.dataset === r.dataset);
    const ok = p && p.swappedSamplesCaught === r.swappedSamplesCaught && p.falseFlagsPer1000 === r.falseFlagsPer1000;
    ok ? same++ : diff++;
    const label = (f.includes('braycurtis') ? '[Bray-Curtis] ' : '') + r.dataset;
    console.log(`${ok ? 'MATCH' : 'DIFF '}  ${label.padEnd(48)} caught ${(100 * r.swappedSamplesCaught).toFixed(1)}%` +
      `  other flags per 1,000: ${r.falseFlagsPer1000}` + (ok ? '' : `  (published: ${p ? (100 * p.swappedSamplesCaught).toFixed(1) + '%, ' + p.falseFlagsPer1000 : 'missing'})`));
  }
}
console.log(`\n${same} of ${same + diff} results match the published numbers exactly.`);
process.exit(diff ? 1 : 0);
JS
