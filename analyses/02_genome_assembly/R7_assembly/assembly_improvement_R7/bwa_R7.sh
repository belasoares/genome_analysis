#!/bin/bash
#SBATCH -A uppmax2024-2-7 -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:30:00
#SBATCH -J alignment_R7
#SBATCH --mem=2G
#SBATCH -o %j.out
#SBATCH -e %j.err

module load bioinfo-tools
module load bwa
module load samtools

READS_DIR=/proj/uppmax2024-2-7/Genome_Analysis/2_Beganovic_2023/DNA_reads/short_reads
ASSEMBLY_DIR=/home/bela/genome_analysis/analyses/02_genome_assembly/R7_assembly/out_nanopore_R7
STORAGE_DIR=/proj/uppmax2024-2-7/nobackup/work/bela_paper2

cd $STORAGE_DIR

bwa index $ASSEMBLY_DIR/assembly.fasta

bwa mem $ASSEMBLY_DIR/assembly.fasta $READS_DIR/SRR24413071_1.fastq.gz $READS_DIR/SRR24413071_2.fastq.gz > aln_R7.sam

samtools view -Sb aln_R7.sam > aln_R7.bam
samtools sort aln_R7.bam -o aln_R7.sorted.bam
rm aln_R7.sam
