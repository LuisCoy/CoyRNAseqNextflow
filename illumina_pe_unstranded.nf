#! /usr/bin/env nextflow

params.reads = "$projectDir/data/reads_{1,2}.fq"
params.outdir = "results" 

log.info """\
    ILLUMINA PR UNSTRANDED WORKFLOW
    -------------------------------
    reads:  ${params.reads}
    outdir: ${params.outdir}
    """
    .stripIndent()

process FASTQC {
    container 'biocontainers/fastqc:v0.11.5'
    
    tag "FastQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc ${reads} -o fastqc_${sample_id}_logs
    """ 
}

workflow {
    Channel 
        .fromFilePairs( params.reads, checkIfExists:true )
        .set { read_pairs_ch }

    fastqc_ch = FASTQC(read_pairs_ch)

}