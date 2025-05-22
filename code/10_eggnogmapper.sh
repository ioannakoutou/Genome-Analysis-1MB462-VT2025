#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -J eggnogmapper_anot
#SBATCH --mail-type=ALL
#SBATCH --mail-user ioanna-kyriaki.koutou.9176@student.uu.se

#modules
module load bioinfo-tools
module load eggNOG-mapper/2.1.9

#
export WORKDIR="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/02_genome_annotation"
export IN="$WORKDIR/braker2_output/braker.agat.aa"
export OUT="$WORKDIR/eggnogmapper_output"
export THREADS=4

#run eggNOGmapper
emapper.py \
  -i "$IN" \
  -o eggnog_anot \
  --itype proteins \
  --output_dir "$OUT" \
  --cpu "$THREADS" \
  --go_evidence experimental \
  --override \
  --decorate_gff "$WORKDIR/braker2_output/braker.cleaned.IDs.gff3" \
  --decorate_gff_ID_field ID
