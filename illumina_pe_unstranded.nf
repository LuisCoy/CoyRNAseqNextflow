#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

params.reads = "$projectDir/data/reads_{1,2}.fq"
params.outdir = "results"
params.genome_indices = "genome_indices"
params.ref_assembly = ""
params.ref_annotations = ""

log.info """\
    ILLUMINA PR UNSTRANDED WORKFLOW
    -------------------------------
    reads:  ${params.reads}
    outdir: ${params.outdir}
    genome indices: ${params.genome_indices} 
    ref assembly = ${params.ref_assembly}
    ref annotations = ${params.ref_annotations}
    """
    .stripIndent()

//Index reference genome
process STAR_GENOME_GEN {
    container 'alexdobin/star:2.6.1d'

    input:
    path ref_assembly
    path ref_annotations

    output:
    path params.genome_indices, emit: genome_index

    script:
    """
    mkdir -p ${params.genome_indices}
    STAR \\
        --runThreadN 14 \\
        --runMode genomeGenerate \\
        --genomeDir ${params.genome_indices} \\
        --genomeFastaFiles ${ref_assembly} \\
        --sjdbGTFfile ${ref_annotations} \\
    -    -sjdbOverhang 100
    """
}
//Fastqc on raw reads
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
//Multiqc on raw reads
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
//Trim raw reads
process TRIMMOMATIC {
    debug true
    container 'quay.io/biocontainers/trimmomatic:0.39--1'
    debug "Using ${task.cpus} CPUs for trimmomatic"
    tag "Trimmomatic on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1P.fq.gz"), path("${sample_id}_2P.fq.gz"), emit: trimmed_reads
    path "trimmomatic_${sample_id}_logs", emit: trimmomatic_logs

    script:
    """
    mkdir trimmomatic_${sample_id}_logs
    trimmomatic \\
        PE \\
        -threads 14 \\
        ${reads} \\
        ${sample_id}_1P.fq.gz \\
        ${sample_id}_1U.fq.gz \\
        ${sample_id}_2P.fq.gz \\
        ${sample_id}_2U.fq.gz \\
        ILLUMINACLIP:TruSeq3-SE:2:30:10 \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:51
    """
}
//Fastqc on trimmed reads
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
//Multiqc on trimmed reads
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
//Star alignment
process STAR_ALIGN{
    container 'alexdobin/star:2.6.1d'
    tag "Star alignment on $sample_id"

    input:
    path genome_indices
    tuple val(sample_id), path(trimmed_read1), path(trimmed_read2)

    output:
    tuple val(sample_id), path("${sample_id}Aligned.sortedByCoord.out.bam"), emit: aligned_sorted
    
    script:
    """
    echo "Debug: Sample ID is $sample_id"
    echo "Debug: Genome indices directory contents:"
    ls -l ${genome_indices}
    echo "Debug: Trimmed read 1 is ${trimmed_read1}"
    echo "Debug: Trimmed read 2 is ${trimmed_read2}"
    echo STAR \\
        --runThreadN 14 \\
        --genomeDir ${genome_indices} \\
        --readFilesIn ${trimmed_read1} ${trimmed_read2} \\
        --readFilesCommand gunzip -c \\
        --outFileNamePrefix ${sample_id}\\
        --outSAMtype BAM SortedByCoordinate

    STAR \\
        --runThreadN 14 \\
        --genomeDir ${genome_indices} \\
        --readFilesIn ${trimmed_read1} ${trimmed_read2} \\
        --readFilesCommand gunzip -c \\
        --outFileNamePrefix ${sample_id}\\
        --outSAMtype BAM SortedByCoordinate
    """
}

//flagstat
process FLAGSTAT {
    tag "Alignment QC on $sample_id"

    publishDir {"${params.outdir}/Alignment_QC"}, mode: "copy"

    input: 
    tuple val(sample_id), path(aligned_bam)

    output:
    path "${sample_id}_flagstat.txt"

    script:
    """
    samtools flagstat ${aligned_bam} > ${sample_id}_flagstat.txt
    """
}

//Count
process HTSEQ_COUNT{

    input:
    tuple val(sample_id), path(aligned_bam)
    path ref_annotations

    output:
    tuple val(sample_id), path("${sample_id}_count.txt"), emit: htseq_count

    script:
    """
    htseq-count \\
        -m union \\
        -s no \\
        -r pos \\
        -a 10 \\
        -i gene_id \\
        -f bam \\
        ${aligned_bam} \\
        ${ref_annotations} \\
        > ${sample_id}_count.txt
    """

}

//Workflow
workflow {
    ref_assembly_ch = channel.fromPath(params.ref_assembly, checkIfExists:true)
    ref_annotations_ch = channel.fromPath(params.ref_annotations,  checkIfExists:true)
    Channel 
        .fromFilePairs( params.reads, checkIfExists:true )
        .set { read_pairs_ch }
    star_genome_gen_ch = STAR_GENOME_GEN(ref_assembly_ch, ref_annotations_ch)
    star_genome_gen_ch.view({"Star genome generator output: $it"})
    
    fastqc_ch = FASTQC(read_pairs_ch)
    MULTIQC(fastqc_ch.collect())
    
    trim_ch = TRIMMOMATIC(read_pairs_ch)
    trim_ch.trimmed_reads.view({"Trimmomatic output: $it"})
    
    fastqc_trimmed_ch = FASTQC_TRIMMED(trim_ch.trimmed_reads)
    MULTIQC_TRIMMED(fastqc_trimmed_ch.collect())
    
    star_align_ch = STAR_ALIGN(star_genome_gen_ch.genome_index, trim_ch.trimmed_reads)
    star_align_ch.view({"Star align output: $it"})

    starqc_ch = FLAGSTAT(star_align_ch.aligned_sorted)
    starqc_ch.view({"Alignment QC output: $it"})

    htseq_count_ch = HTSEQ_COUNT(star_align_ch.aligned_sorted, ref_annotations_ch)
    htseq_count_ch.view({"htseq count output: $it"})
}