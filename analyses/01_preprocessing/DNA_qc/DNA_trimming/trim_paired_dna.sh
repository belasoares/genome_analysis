#!/bin/bash
#SBATCH -A uppmax2024-2-7 -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:45:00
#SBATCH -J trim_paired_dna
#SBATCH --mem=6G
#SBATCH -o paireddna_%j.out
#SBATCH -e paireddna_%j.err

module load bioinfo-tools
module load trimmomatic

INPUT_DIR=/proj/uppmax2024-2-7/Genome_Analysis/2_Beganovic_2023/DNA_reads/short_reads
OUTPUT_DIR=/home/bela/genome_analysis/analyses/01_preprocessing/trimming

for f1 in $INPUT_DIR/*_1.fastq.gz
do
    f2=${f1%%_1.fastq.gz}"_2.fastq.gz" 
    base_f1=$(basename $f1) 
    base_f2=$(basename $f2) 
    
    trimmomatic PE $f1 $f2 $OUTPUT_DIR/forward_paired_${base_f1} $OUTPUT_DIR/forward_unpaired_${base_f1} $OUTPUT_DIR/reverse_paired_${base_f2} $OUTPUT_DIR/reverse_unpaired_${base_f2} LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:100
done



