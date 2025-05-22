#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:40:00
#SBATCH -J FastQC
#SBATCH --mail-type=ALL
#SBATCH --mail-user ioanna-kyriaki.koutou.9176@student.uu.se

#modules
module load bioinfo-tools
module load FastQC/0.11.9

#
export IN1="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/data/dna_seq/illumina/raw/chr3_illumina_R1.fastq.gz"
export IN2="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/data/dna_seq/illumina/raw/chr3_illumina_R2.fastq.gz"
export OUT="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/00_read_preprocessing/fastqc_output"
export THREADS=2

#run FastQC
fastqc -t "$THREADS" -o "$OUT" "$IN1" "$IN2"

