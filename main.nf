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

process arcasHLA_extract {
    container "${params.container__arcashla}"
    publishDir "${params.outdir}/${sample}", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("hla/*.extracted.1.fq.gz"), path("hla/*.extracted.2.fq.gz"), emit: fq
    tuple val(sample), path("*.extract.log"), emit: log

    """#!/bin/bash
set -e

arcasHLA extract \
    --threads ${task.cpus} \
    --outdir hla \
    --verbose \
    --temp tmp \
    "${bam}" 2>&1 | tee "${sample}.extract.log"

ls -lahtr

for NUM in 1 2; do
    if [ ! -s "hla/${bam.name.replaceAll(/.bam/, '')}.extracted.\$NUM.fq.gz" ]; then
        gzip "hla/${bam.name.replaceAll(/.bam/, '')}.extracted.\$NUM.fq"
    fi
done

"""
}

process arcasHLA_genotype {
    container "${params.container__arcashla}"
    publishDir "${params.outdir}/${sample}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(fq1), path(fq2), path(db_tar)

    output:
    path "hla/*"

    """#!/bin/bash
set -e

# Unpack the database tar file
unpack_database.sh ${db_tar}

echo "\$(date) - Running genotype"
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
    "${fq1}" "${fq2}" 2>&1 | tee "${sample}.genotype.log"
"""
}

process merge_arcasHLA_results {
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

    // Get a tar of the reference database
    imgthla_tar = Channel.fromPath(
        "${params.database_dir}/v${params.database}-alpha.tar.gz",
        checkIfExists: true
    )

    arcasHLA_extract(STAR.out.mix(bam_input))

    arcasHLA_genotype(arcasHLA_extract.out.fq.combine(imgthla_tar))

    merge_arcasHLA_results(arcasHLA_genotype.out.flatten().toSortedList())
}