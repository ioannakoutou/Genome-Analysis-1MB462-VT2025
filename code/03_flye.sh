#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH -J Flye_Assembly
#SBATCH --mail-type=ALL
#SBATCH --mail-user ioanna-kyriaki.koutou.9176@student.uu.se

#modules
module load bioinfo-tools
module load Flye/2.9.5
 
#
export IN="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/data/dna_seq/nanopore/chr3_clean_nanopore.fq.gz"
export OUT="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/01_genome_assembly/flye_output"
export THREADS=8


# run flye
flye --nano-raw "$IN" --threads "$THREADS" --out-dir "$OUT" --resume
