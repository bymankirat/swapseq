# Host-DNA check (SwapSeq 1.2)

Shotgun reads from stool, swabs or saliva always carry some of the person's own DNA. Genotypes at common
human SNPs identify a person, however much their microbes change. That covers the cases where the microbial
fingerprint fails: newborns, heavy antibiotic use, and acute illness.

SwapSeq scores each sample's host reads against every person's pooled reads from all their other samples,
the same idea that single-cell demultiplexing uses. The person whose genotypes best explain the reads wins.
A sample is flagged by host DNA when another person explains its reads at least e^8 (about 3,000) times
better than the person on its label.

**Status.** Tested on simulated genotypes only, not yet on real samples. With a panel of 200,000 common SNPs
and one read per covered site (40 people, 6 samples each), it caught 95% of swaps at 2,000 covered sites per
sample and all of them from 2,500, with no false flags; at 1,500 sites it caught 25% and made one false flag in
about 1,200 samples. SwapSeq therefore judges only samples with at least 2,000 covered sites. Treat results on
real data as a strong lead and confirm them against lab records.

## Before you start: privacy and ethics

Host reads and genotypes are **human genetic data**, a special category of personal data under the GDPR and
similar laws. Only use this check on data your ethics approval and consent cover. Keep the files on your own
systems: SwapSeq runs locally and never uploads anything. Do not share or publish the host count file.
Public microbiome data usually has host reads removed, so this check is for a lab's own raw data.

## 1. A panel of common SNPs (once)

Use common, biallelic SNPs, for example from the 1000 Genomes sites on GRCh38:

```bash
bcftools view -m2 -M2 -v snps -i 'INFO/AF>=0.05 && INFO/AF<=0.95' 1kGP_sites.GRCh38.vcf.gz -Oz -o panel.vcf.gz
bcftools index panel.vcf.gz
```

More sites mean more overlap between samples. A few million common SNPs work better than a small panel.

## 2. Host reads per sample

Most shotgun pipelines already map reads to the human genome to remove them, for example with KneadData
or `bowtie2 -x GRCh38`. Keep the reads that **did** map (the host BAM) instead of throwing them away:

```bash
bowtie2 -x GRCh38 -1 S1_R1.fq.gz -2 S1_R2.fq.gz --no-unal | samtools sort -o S1.host.bam
samtools index S1.host.bam
```

## 3. Allele counts at the panel sites

```bash
echo -e "sample\tsite\talt\ttotal" > host_counts.tsv
for bam in *.host.bam; do
  s=$(basename "$bam" .host.bam)
  bcftools mpileup -f GRCh38.fa -R panel.vcf.gz -a AD -d 50 -Ou "$bam" \
    | bcftools query -f '%CHROM:%POS[\t%AD]\n' \
    | awk -v s="$s" -F'\t' '{split($2,a,","); t=a[1]+a[2]; if (t>0) print s"\t"$1"\t"a[2]"\t"t}' >> host_counts.tsv
done
```

The sample names must match the sample IDs in your feature table.

## 4. Run SwapSeq with the host counts

Command line:

```bash
npx github:bymankirat/swapseq table.tsv metadata.tsv --subject subject_id --host host_counts.tsv
```

R:

```r
host <- read.delim("host_counts.tsv")
res <- swapseqr::check_swaps(ps, subject = "subject_id", host = host)
res[res$flagged_by_host, c("sample", "subject", "host_looks_like", "host_llr", "host_sites")]
```

The columns `host_verdict` (same person, other person, unclear or unknown), `host_looks_like`, `host_llr`
(log likelihood ratio; above 8 means another person fits much better) and `host_sites` show the host
evidence for every sample. A sample is judged only when it has at least 2,000 covered sites.
