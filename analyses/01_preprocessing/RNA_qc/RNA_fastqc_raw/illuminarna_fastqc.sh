#!/bin/bash
#SBATCH -A uppmax2024-2-11 -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:30:00
#SBATCH -J illuminarna_fastqc
#SBATCH --mem=6G
#SBATCH -o %j.out
#SBATCH -e %j.err

module load bioinfo-tools
module load FastQC/0.11.9

INPUT_DIR=/proj/uppmax2024-2-7/Genome_Analysis/2_Beganovic_2023/RNA_reads
OUTPUT_DIR=/home/bela/genome_analysis/analyses/01_preprocessing/RNA_qc/RNA_fastqc_raw

fastqc -o $OUTPUT_DIR --noextract -t 6 \
$INPUT_DIR/SRR24516459_1.fastq.gz \
$INPUT_DIR/SRR24516459_2.fastq.gz \
$INPUT_DIR/SRR24516460_1.fastq.gz \
$INPUT_DIR/SRR24516460_2.fastq.gz \
$INPUT_DIR/SRR24516461_1.fastq.gz \
$INPUT_DIR/SRR24516461_2.fastq.gz
