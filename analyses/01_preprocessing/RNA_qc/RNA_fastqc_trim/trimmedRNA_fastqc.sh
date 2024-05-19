#!/bin/bash
#SBATCH -A uppmax2024-2-7 -M snowy
#SBATCH -p core                   
#SBATCH -n 6
#SBATCH -t 00:45:00                
#SBATCH -J trimmedRNA_fastqc         
#SBATCH --mem=6G                     
#SBATCH -o %j.out       
#SBATCH -e %j.err

module load bioinfo-tools
module load FastQC/0.11.9

INPUT_DIR=/home/bela/genome_analysis/analyses/01_preprocessing/RNA_qc/RNA_trimming
OUTPUT_DIR=/home/bela/genome_analysis/analyses/01_preprocessing/RNA_qc/RNA_fastqc_trim

fastqc -o $OUTPUT_DIR --noextract -t 6 $INPUT_DIR/*.gz
