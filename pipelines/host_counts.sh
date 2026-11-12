#!/usr/bin/env bash
# host_counts.sh - make host_counts.tsv for SwapSeq's host-DNA check from the human reads of each sample.
#
#   pipelines/host_counts.sh -b BAM_DIR -p panel.vcf.gz -f GRCh38.fa [-o host_counts.tsv] [-s .host.bam] [-d 50]
#
#   -b  folder with one BAM per sample: the reads that mapped to the human genome (sorted)
#   -p  panel of common SNPs (bgzipped VCF with index), see docs/host-dna.md
#   -f  the human reference the BAMs were mapped to (FASTA)
#   -o  output file (default host_counts.tsv)
#   -s  BAM file ending; the sample ID is the file name without it (default .host.bam)
#   -d  maximum depth per site (default 50)
#
# Needs samtools and bcftools. Writes OUTPUT (sample, site, alt, total) and OUTPUT.summary.tsv (covered sites per
# sample). Host reads are human genetic data: use only what your ethics approval and consent cover, keep the files on
# your own systems and do not share them. Sample IDs must match the IDs in your feature table.
set -euo pipefail
out=host_counts.tsv; suffix=.host.bam; depth=50; bamdir=""; panel=""; ref=""
usage(){ sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-1}"; }
while getopts "b:p:f:o:s:d:h" o; do case "$o" in
  b) bamdir=$OPTARG;; p) panel=$OPTARG;; f) ref=$OPTARG;; o) out=$OPTARG;; s) suffix=$OPTARG;; d) depth=$OPTARG;; h) usage 0;; *) usage 1;; esac; done
die(){ echo "host_counts: $*" >&2; exit 1; }
[ -n "$bamdir" ] && [ -n "$panel" ] && [ -n "$ref" ] || usage 1
for t in samtools bcftools; do command -v "$t" >/dev/null 2>&1 || die "$t is not installed (needed to read the BAM files)"; done
[ -d "$bamdir" ] || die "no such folder: $bamdir"
[ -f "$panel" ] || die "panel not found: $panel"
[ -f "$panel.tbi" ] || [ -f "$panel.csi" ] || die "the panel has no index; run: bcftools index $panel"
[ -f "$ref" ] || die "reference not found: $ref"
case "$depth" in ''|*[!0-9]*) die "-d must be a whole number";; esac
shopt -s nullglob; bams=("$bamdir"/*"$suffix"); shopt -u nullglob
[ "${#bams[@]}" -gt 0 ] || die "no files ending in $suffix in $bamdir (use -s to give another ending)"
for b in "${bams[@]}"; do [ "$b" -ef "$out" ] && die "the output would overwrite an input file"; done
tmp=$(mktemp "${out}.XXXXXX"); sum="${out}.summary.tsv"; trap 'rm -f "$tmp"' EXIT
printf 'sample\tsite\talt\ttotal\n' > "$tmp"; printf 'sample\tcovered_sites\tenough\n' > "$sum"
low=0; n=0
for bam in "${bams[@]}"; do
  s=$(basename "$bam" "$suffix"); n=$((n+1))
  [ -f "$bam.bai" ] || [ -f "$bam.csi" ] || [ -f "${bam%.bam}.bai" ] || samtools index "$bam"
  sites=$(bcftools mpileup -f "$ref" -R "$panel" -a AD -d "$depth" -Ou "$bam" 2>/dev/null \
    | bcftools query -f '%CHROM:%POS[\t%AD]\n' \
    | awk -v s="$s" -F'\t' 'BEGIN{c=0} {k=split($2,a,","); if(k<2) next; t=a[1]+a[2]; if(t>0){print s"\t"$1"\t"a[2]"\t"t; c++}} END{print c > "/dev/stderr"}' 2>&1 >> "$tmp" | tail -n 1)
  sites=${sites:-0}; ok=yes; if [ "$sites" -lt 1500 ]; then ok=no; low=$((low+1)); fi
  printf '%s\t%s\t%s\n' "$s" "$sites" "$ok" >> "$sum"
  echo "  $s: $sites covered sites$([ "$ok" = no ] && echo '  (under 1,500: will not be judged)')" >&2
done
mv "$tmp" "$out"; trap - EXIT
echo "host_counts: wrote $out for $n samples; $low with fewer than 1,500 covered sites. Summary: $sum" >&2
echo "Next: npx github:bymankirat/swapseq table.tsv metadata.tsv --subject SUBJECT_COLUMN --host $out" >&2
