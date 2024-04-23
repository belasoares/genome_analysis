#!/bin/bash
#SBATCH -A uppmax2024-2-7 -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:30:00
#SBATCH -J pilon_HP126
#SBATCH -o %j.out
#SBATCH -e %j.err

module load bioinfo-tools
module load samtools
module load Pilon/1.24

ASSEMBLY_DIR=/home/bela/genome_analysis/analyses/02_genome_assembly/HP126_assembly/out_nanopore_HP126
STORAGE_DIR=/proj/uppmax2024-2-7/nobackup/work/bela_paper2

samtools index $STORAGE_DIR/aln_HP126.sorted.bam

java -jar $PILON_HOME/pilon.jar --genome $ASSEMBLY_DIR/assembly.fasta --bam $STORAGE_DIR/aln_HP126.sorted.bam --output HP126
