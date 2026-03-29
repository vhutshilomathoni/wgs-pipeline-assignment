#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// 1. Quality contrl
process FASTQC {
    tag "QC on ${id}"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(id), path(reads)

    output:
    path "*.{html,zip}"

    script:
    """
    fastqc ${reads}
    """
}

// 2. trimming (using my singularity container)
process TRIMMOMATIC {
    tag "Trimming ${id}"
    publishDir "${params.outdir}/trim", mode: 'copy'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_R{1,2}_paired.fq.gz"), emit: trimmed_reads

    script:
    """
    java -jar /opt/trimmomatic/trimmomatic-0.39.jar PE \
    ${reads[0]} ${reads[1]} \
    ${id}_R1_paired.fq.gz ${id}_R1_unpaired.fq.gz \
    ${id}_R2_paired.fq.gz ${id}_R2_unpaired.fq.gz \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// 3. Alignment
process BWA {
    publishDir "${params.outdir}/alignment", mode: 'copy'
    tag "Aligning ${id}"
    
    input:
    path genome_fasta
    path genome_indices
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}.bam")

    script:
    """
    bwa mem -t ${task.cpus} ${genome_fasta} ${reads} | samtools view -Sb - > ${id}.bam
    """
}

// 4. variant calling
process VARIANT_CALLING {
    tag "Calling variants for ${id}"
    publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    path genome_fasta
    path genome_indices
    tuple val(id), path(bam)

    output:
    path "${id}.vcf", emit: vcf

    script:
    """
    samtools sort ${bam} -o sorted.bam
    samtools index sorted.bam
    bcftools mpileup -f ${genome_fasta} sorted.bam | bcftools call -mv -Ov -o ${id}.vcf
    """
}

process CREATE_DB {
    publishDir "${params.outdir}/db", mode: 'copy'

    input:
    path vcf

    output:
    path "variants.db"

    script:
    """
    #1. extrcting columns 1, 2, 4, 5, AND 6
    grep -v "^#" ${vcf} | awk 'BEGIN {OFS=","} {print \$1, \$2, \$4, \$5, \$6}' > variants.csv

   
    sqlite3 variants.db <<EOF
.mode csv
CREATE TABLE variants (chrom TEXT, pos INT, ref TEXT, alt TEXT, qual REAL);
.import variants.csv variants
EOF
    """
}

// Main Workflow Logic
workflow {
    // Create input channels
    read_pairs_ch = Channel.fromFilePairs(params.fastq, checkIfExists: true)
    genome_fasta = file(params.genome)
    genome_indices = Channel.fromPath("${params.genome}.*").collect()

    // Run processes
    FASTQC(read_pairs_ch)
    TRIMMOMATIC(read_pairs_ch)
    
    // Pass genome AND indices to alignment/calling
    BWA(genome_fasta, genome_indices, TRIMMOMATIC.out.trimmed_reads)
    VARIANT_CALLING(genome_fasta, genome_indices, BWA.out)
    
    // Create the final SQLite database
    CREATE_DB(VARIANT_CALLING.out.vcf)
}
