#!/bin/bash
#SBATCH -A uppmax2024-2-7 -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:10:00
#SBATCH -J quast_HP126
#SBATCH -o %j.out
#SBATCH -e %j.err

module load bioinfo-tools
module load quast

ASSEMBLY_DIR=/home/bela/genome_analysis/analyses/02_genome_assembly/HP126_assembly/assembly_improvement_HP126
REF_DIR=/proj/uppmax2024-2-7/Genome_Analysis/2_Beganovic_2023/reference_genome

quast.py -r $REF_DIR/HP126_genome.fasta --large $ASSEMBLY_DIR/HP126.fasta 
