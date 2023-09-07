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

STAR --genomeDir ${genomeDir} \
     --readFilesIn ${R1} ${R2} \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --runThreadN ${task.cpus} \
     --outTmpDir tmp \
     --outFileNamePrefix "${sample}."

    """
}

process arcasHLA {
    container "${params.container__arcashla}"
    publishDir "${params.outdir}/${sample}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(bam)

    output:
    path "hla/*.genotype.json", emit: json
    path "hla/"

    """#!/bin/bash
set -e

arcasHLA extract \
    --threads ${task.cpus} \
    --outdir hla \
    --verbose \
    --temp tmp \
    "${bam}"

ls -lahtr

for NUM in 1 2; do
    if [ ! -s "hla/${bam.replaceAll(/.bam/, '')}.extracted.\$NUM.fq.gz" ]; then
        gzip "hla/${bam.replaceAll(/.bam/, '')}.extracted.\$NUM.fq"
    fi
done

# Required to install the database
git lfs install

# Set up the reference database
arcasHLA reference --version "${params.database}"

arcasHLA genotype \
    --genes ${params.genes} \
    --population "${params.population}" \
    --zygosity_threshold "${params.zygosity_threshold}" \
    --min_count "${params.min_count}" \
    --outdir hla \
    --threads ${task.cpus} \
    --verbose \
    --temp tmp \
    --log "hla/${sample}.log" \
    "hla/${bam.replaceAll(/.bam/, '')}.extracted.1.fq.gz" "hla/${bam.replaceAll(/.bam/, '')}.extracted.2.fq.gz"
"""
}

process merge {
    container "${params.container__arcashla}"
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true

    input:
    path "*"

    output:
    path "*"

    """#!/bin/bash
set -e

arcasHLA merge
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
        .branch {
            fastq: it.fastq_1 != null
            bam: true
        }
        .set { inputs }

    inputs
        .fastq
        .map {row -> [
            row.sample,
            file(row.fastq_1, checkIfExists: true),
            file(row.fastq_2, checkIfExists: true)
        ]}
        .set { fastq_input }

    inputs
        .bam
        .map {row -> [
            row.sample,
            file(row.bam)
        ]}
        .set { bam_input }

    genomeDir = file(params.genomeDir, checkIfExists: true, type: 'dir')

    STAR(fastq_input, genomeDir)

    arcasHLA(STAR.out.mix(bam_input))

    merge(arcasHLA.out.json.toSortedList())
}