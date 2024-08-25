#!/bin/bash
#SBATCH -A uppmax2024-2-11 -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:45:00
#SBATCH -J trim_paired_RNA
#SBATCH -o %j.out
#SBATCH -e %j.err

module load bioinfo-tools
module load trimmomatic

INPUT_DIR=/proj/uppmax2024-2-7/Genome_Analysis/2_Beganovic_2023/RNA_reads
OUTPUT_DIR=/home/bela/genome_analysis/analyses/01_preprocessing/RNA_qc/RNA_trimming

declare -a samples=("SRR24516459" "SRR24516462" "SRR24516463" "SRR24516464")

for sample_id in "${samples[@]}"
do
    f1=$INPUT_DIR/${sample_id}_1.fastq.gz
    f2=$INPUT_DIR/${sample_id}_2.fastq.gz
    base_f1=$(basename $f1)
    base_f2=$(basename $f2)

    trimmomatic PE $f1 $f2 \
    $OUTPUT_DIR/forward_paired_${base_f1} \
    $OUTPUT_DIR/forward_unpaired_${base_f1} \
    $OUTPUT_DIR/reverse_paired_${base_f2} \
    $OUTPUT_DIR/reverse_unpaired_${base_f2} \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:20
done
