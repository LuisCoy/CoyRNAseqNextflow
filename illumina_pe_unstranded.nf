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
    cpus 8
    container 'quay.io/biocontainers/trimmomatic:0.39--1'
    debug "Using ${task.cpus} CPUs for trimmomatic"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1P.fq.gz"), path("${sample_id}_2P.fq.gz"), emit: trimmed_reads
    path "trimmomatic_${sample_id}_logs"

    script:
    """
    mkdir trimmomatic_${sample_id}_logs
    trimmomatic \\
        PE \\
        -threads ${task.cpus} \\
        ${reads} \\
        ${sample_id}_1P.fq.gz \\
        ${sample_id}_1U.fq.gz \\
        ${sample_id}_2P.fq.gz \\
        ${sample_id}_2U.fq.gz \\
        ILLUMINACLIP:TruSeq3-SE:2:30:10 \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:51
    """
}

process FASTQC_TRIMMED {
    container 'biocontainers/fastqc:v0.11.5'
    publishDir {"${params.outdir}/trimmed"}, mode:'copy'
    tag "FastQC on $sample_id"

    input:
    tuple val(sample_id), path(trimmed_read1), path(trimmed_read2)

    output:
    path "fastqc_trimmed_${sample_id}_logs"

    script:
    """
    mkdir fastqc_trimmed_${sample_id}_logs
    fastqc ${trimmed_read1} ${trimmed_read2} \\
    -o fastqc_trimmed_${sample_id}_logs \\
    --threads ${task.cpus}
    """ 
}

process MULTIQC_TRIMMED {
    container 'multiqc/multiqc:latest'

    publishDir {"${params.outdir}/trimmed"}, mode:'copy'
    
    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow {
    Channel 
        .fromFilePairs( params.reads, checkIfExists:true )
        .set { read_pairs_ch }

    fastqc_ch = FASTQC(read_pairs_ch)
    MULTIQC(fastqc_ch.collect())
    trim_ch = TRIMMOMATIC(read_pairs_ch)
    fastqc_trimmed_ch = FASTQC_TRIMMED(trim_ch.trimmed_reads)
    MULTIQC_TRIMMED(fastqc_trimmed_ch.collect())
}