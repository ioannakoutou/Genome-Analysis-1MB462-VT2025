#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 4-00:00:00
#SBATCH -J Canu_Assembly
#SBATCH --mail-type=ALL
#SBATCH --mail-user ioanna-kyriaki.koutou.9176@student.uu.se

#modules
module load bioinfo-tools
module load canu/2.2

#
IN="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/data/dna_seq/nanopore/chr3_clean_nanopore.fq.gz"
OUT="/proj/uppmax2025-3-3/nobackup/work/ioko9176/canu_output"
PREFIX="canu_chr3"
GENOMESIZE="17m"

# run canu
canu -p "$PREFIX" -d "$OUT" genomeSize="$GENOMESIZE" -nanopore "$IN" useGrid=false
