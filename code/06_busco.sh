#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 00:30:00
#SBATCH -J busco_assessment
#SBATCH --mail-type=ALL
#SBATCH --mail-user ioanna-kyriaki.koutou.9176@student.uu.se

#modules
module load bioinfo-tools
module load BUSCO/5.7.1

#
export WORKDIR="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/01_genome_assembly"
export IN="$WORKDIR/pilon_output/polished_assembly.fasta"
export LINEAGE="embryophyta_odb10"
export THREADS=4

cd "$WORKDIR"

# download lineage dataset
busco --download "$LINEAGE" --download_path "$WORKDIR/busco_downloads"

#run busco
busco -i "$IN" -o "busco_output" -l "$WORKDIR/busco_downloads/lineages/$LINEAGE" -m genome -c $THREADS --out_path "$WORKDIR" -f
