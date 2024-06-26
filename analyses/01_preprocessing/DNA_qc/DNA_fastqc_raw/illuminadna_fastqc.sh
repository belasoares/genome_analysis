#!/bin/bash
#SBATCH -A uppmax2024-2-7 -M snowy
#SBATCH -p core                   
#SBATCH -n 6
#SBATCH -t 00:45:00                
#SBATCH -J illuminadna_fastqc         
#SBATCH --mem=6G                     
#SBATCH -o %j.out       
#SBATCH -e %j.err

module load bioinfo-tools
module load FastQC/0.11.9

INPUT_DIR=/proj/uppmax2024-2-7/Genome_Analysis/2_Beganovic_2023/DNA_reads/short_reads
OUTPUT_DIR=/home/bela/genome_analysis/analyses/01_preprocessing/fastqc_raw

fastqc -o $OUTPUT_DIR --noextract -t 6 $INPUT_DIR/*.gz


