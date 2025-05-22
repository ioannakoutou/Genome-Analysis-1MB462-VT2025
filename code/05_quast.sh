#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J quast_assessment
#SBATCH --mail-type=ALL
#SBATCH --mail-user ioanna-kyriaki.koutou.9176@student.uu.se

#modules
module load bioinfo-tools
module load quast/5.0.2

#
export WORKDIR="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/01_genome_assembly"
export IN="$WORKDIR/pilon_output/polished_assembly.fasta"
export OUT="$WORKDIR/quast_output"
export THREADS=2

#run quast
quast.py "$IN" -o "$OUT" --threads "$THREADS"

