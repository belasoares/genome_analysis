#!/bin/bash
#SBATCH -A uppmax2024-2-7 -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 04:00:00
#SBATCH -J HTSeq_count
#SBATCH -o %j.out
#SBATCH -e %j.err

module load bioinfo-tools
module load python
module load samtools/1.20

FEATURES_FILE="/home/bela/genome_analysis/analyses/03_structural_annotation/HP126/HP126without.gff"
INPUT_DIR="/home/bela/genome_analysis/analyses/04_RNA_mapping"
OUTPUT_DIR="/home/bela/genome_analysis/analyses/05_expression_analysis/read_count"
OUTPUT_FILE="${OUTPUT_DIR}/gene_counts.txt"

for bam_file in ${INPUT_DIR}/*.bam; do

    base_name=$(basename ${bam_file} .bam)
    
    htseq-count -f bam -r pos -s no -i ID -t CDS ${bam_file} $FEATURES_FILE > ${OUTPUT_DIR}/${base_name}_counts.txt

done
