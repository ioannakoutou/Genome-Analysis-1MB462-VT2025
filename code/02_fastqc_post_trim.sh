#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:30:00
#SBATCH -J FastQC_trimmed
#SBATCH --mail-type=ALL
#SBATCH --mail-user ioanna-kyriaki.koutou.9176@student.uu.se

#modules
module load bioinfo-tools
module load FastQC/0.11.9

#
export IN1="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/00_read_preprocessing/trimmomatic_PE_output/chr3_illumina_R1_trimmed_paired.fastq.gz"
export IN2="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/00_read_preprocessing/trimmomatic_PE_output/chr3_illumina_R2_trimmed_paired.fastq.gz"
export OUT="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/00_read_preprocessing/fastqc_trimmed_output"
export THREADS=2

#run FastQC
fastqc -t "$THREADS" -o "$OUT" "$IN1" "$IN2"
