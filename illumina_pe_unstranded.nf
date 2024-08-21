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

    publishDir params.outdir, mode:'copy'

    debug "Using ${task.cpus} CPUs for fastqc"

    tag "FastQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc ${reads} -o fastqc_${sample_id}_logs --threads ${task.cpus}
    """ 
}

process MULTIQC {
    container 'multiqc/multiqc:latest'
    
    publishDir params.outdir, mode:'copy'
    
    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

process TRIMMOMATIC {
    container 'quay.io/biocontainers/trimmomatic:0.39--1'

    debug "Using ${task.cpus} CPUs for trimmomatic"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "trimmomatic_${sample_id}_logs"

    script:
    """
    mkdir trimmomatic_${sample_id}_logs
    trimmomatic PE -threads ${task.cpus} ${reads} trimmomatic_${sample_id}_logs/${sample_id}_1P.fq.gz trimmomatic_${sample_id}_logs/${sample_id}_1U.fq.gz trimmomatic_${sample_id}_logs/${sample_id}_2P.fq.gz trimmomatic_${sample_id}_logs/${sample_id}_2U.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:51
    """
}

workflow {
    Channel 
        .fromFilePairs( params.reads, checkIfExists:true )
        .set { read_pairs_ch }

    fastqc_ch = FASTQC(read_pairs_ch)
    MULTIQC(fastqc_ch.collect())
    trim_ch = TRIMMOMATIC(read_pairs_ch)
}