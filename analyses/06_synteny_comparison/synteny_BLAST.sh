#!/bin/bash
#SBATCH -A uppmax2024-2-7 -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:45:00
#SBATCH -J synteny_BLAST
#SBATCH -o %j.out
#SBATCH -e %j.err

module load bioinfo-tools
module load blast/2.15.0+
module load java

REF_GENOME="/home/bela/genome_analysis/raw_data/reference_genome/R7_genome.fasta"
QUERY_SEQ="/home/bela/genome_analysis/analyses/02_genome_assembly/HP126_assembly/assembly_improvement_HP126/HP126.fasta"
BLAST_DB="custom_db"
BLAST_OUTPUT="blast_results.txt"

makeblastdb -in $REF_GENOME -dbtype nucl -out $BLAST_DB

blastn -query $QUERY_SEQ -db $BLAST_DB -outfmt 6 -out $BLAST_OUTPUT -evalue 1e-10





