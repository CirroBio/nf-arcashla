#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process STAR {
    container "${params.container__star}"
    publishDir "${params.outdir}/${sample}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(R1), path(R2)
    path genomeDir

    output:
    tuple val(sample), path("*Aligned.sortedByCoord.out.bam")
    
    """#!/bin/bash
set -e

mkdir tmp

STAR --genomeDir ${genomeDir} \
     --readFilesIn ${R1} ${R2} \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --runThreadN ${task.cpus} \
     --outTmpDir tmp \
     --outFileNamePrefix "${sample}"

rm -r tmp
    
    """
}

process arcasHLA {
    container "${params.container__arcashla}"
    publishDir "${params.outdir}/${sample}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(bam)

    output:
    path "hla/"

    """#!/bin/bash
set -e

mkdir tmp

arcasHLA extract \
    --threads ${task.cpus} \
    --paired \
    --outdir hla \
    --verbose \
    --temp tmp \
    "${bam}"

arcasHLA genotype \
    --genes ${params.genes} \
    --outdir hla \
    --threads ${task.cpus} \
    --verbose \
    --temp tmp \
    "hla/${sample}.Aligned.sortedByCoord.out.extracted.1.fq.gz" "hla/${sample}.Aligned.sortedByCoord.out.extracted.2.fq.gz"

rm -r tmp
"""
}

workflow {
    // Parse the input files from the samplesheet
    Channel
        .fromPath(
            params.samplesheet,
            checkIfExists: true
        )
        .splitCsv(header: true)
        .flatten()
        .map {row -> [
            row.sample,
            file(row.fastq_1, checkIfExists: true),
            file(row.fastq_2, checkIfExists: true)
        ]}
        .set { inputs }

    genomeDir = file(params.genomeDir, checkIfExists: true, type: 'dir')

    STAR(inputs, genomeDir)

    arcasHLA(STAR.out)
}