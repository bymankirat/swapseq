// SwapSeq as a Nextflow process: checks every sample against the person on its label before analysis.
// Needs Node.js 18+ on the compute node (or a container with Node). Example:
//   nextflow run main.nf --table feature-table.tsv --metadata metadata.tsv --subject subject_id
// The run stops (exit code 2) when a sample is strongly flagged; drop --fail-on-flags to only report.

params.table    = null
params.metadata = null
params.subject  = 'subject_id'
params.type     = 'counts'          // or 'markers' for strain-marker tables
params.outdir   = 'swapseq'

process SWAPSEQ {
    tag "${table.baseName}"
    publishDir params.outdir, mode: 'copy'

    input:
    path table
    path metadata

    output:
    path 'swapseq.tsv',  emit: audit
    path 'swapseq.json', emit: summary

    script:
    """
    npx --yes github:bymankirat/swapseq ${table} ${metadata} \\
        --subject ${params.subject} --type ${params.type} --out swapseq --fail-on-flags
    """
}

workflow {
    SWAPSEQ(file(params.table, checkIfExists: true), file(params.metadata, checkIfExists: true))
}
