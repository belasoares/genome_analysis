#!/bin/bash
#SBATCH -A uppmax2024-2-7 -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 04:00:00
#SBATCH -J BWA_mapping
#SBATCH -o %j.out
#SBATCH -e %j.err

module load bioinfo-tools
module load bwa/0.7.18
module load samtools

OUTPUT_DIR=/home/bela/genome_analysis/analyses/04_RNA_mapping
REF_GENOME=/home/bela/genome_analysis/analyses/02_genome_assembly/HP126_assembly/assembly_improvement_HP126/HP126.fasta

bwa index $REF_GENOME

READ_DIR=/home/bela/genome_analysis/raw_data/RNA_reads

for f1 in $READ_DIR/*_1.fastq.gz; do
    f2="${f1%_1.fastq.gz}_2.fastq.gz"

    base_name=$(basename "${f1%_1.fastq.gz}")

    output_sam="$OUTPUT_DIR/${base_name}_mapped_reads.sam"
    output_bam="$OUTPUT_DIR/${base_name}_sorted.bam"
   
    bwa mem -t 2 $REF_GENOME "$f1" "$f2" > "$output_sam"

    samtools view -bS "$output_sam" | samtools sort -o "$output_bam"
    samtools index "$output_bam"

    rm "$output_sam"
done
