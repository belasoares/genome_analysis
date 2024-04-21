#!/bin/bash
#SBATCH -A uppmax2024-2-7 -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 02:00:00
#SBATCH -J nanopore_assembly
#SBATCH -o %j.out
#SBATCH -e %j.err

module load bioinfo-tools
module load Flye/2.9.1

INPUT_DIR=/proj/uppmax2024-2-7/Genome_Analysis/2_Beganovic_2023/DNA_reads/SRR24413066.fastq.gz

flye --nano-raw $INPUT_DIR --out-dir out_nanopore_HP126 --threads 4
