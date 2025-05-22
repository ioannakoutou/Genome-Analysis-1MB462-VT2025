#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:40:00
#SBATCH -J Trimmomatic_PE
#SBATCH --mail-type=ALL
#SBATCH --mail-user ioanna-kyriaki.koutou.9176@student.uu.se

#modules
module load bioinfo-tools
module load trimmomatic/0.39

#
export IN1="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/data/dna_seq/illumina/raw/chr3_illumina_R1.fastq.gz"
export IN2="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/data/dna_seq/illumina/raw/chr3_illumina_R2.fastq.gz"
export OUT="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/00_read_preprocessing/trimmomatic_PE_output"
export THREADS=2

#output filenames
export OUT_P1="$OUT/chr3_illumina_R1_trimmed_paired.fastq.gz"
export OUT_UP1="$OUT/chr3_illumina_R1_trimmed_unpaired.fastq.gz"
export OUT_P2="$OUT/chr3_illumina_R2_trimmed_paired.fastq.gz"
export OUT_UP2="$OUT/chr3_illumina_R2_trimmed_unpaired.fastq.gz"

#run trimmomatic PE
trimmomatic PE -threads "$THREADS" -phred33 -trimlog "$OUT/trimmomatic_PE.log" \
  "$IN1" "$IN2" \
  "$OUT_P1" "$OUT_UP1" "$OUT_P2" "$OUT_UP2" \
  ILLUMINACLIP:/sw/bioinfo/trimmomatic/0.39/rackham/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40
