#!/bin/bash
#SBATCH -A uppmax2024-2-7 -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:30:00
#SBATCH -J prokka_HP126.sh
#SBATCH -o %j.out
#SBATCH -e %j.err

ASSEMBLY_DIR=/home/bela/genome_analysis/analyses/02_genome_assembly/HP126_assembly/assembly_improvement_HP126/HP126.fasta

module load bioinfo-tools
module load prokka

prokka --prefix HP126 --genus Streptomyces --species rimosus --strain HP126 --usegenus $ASSEMBLY_DIR
