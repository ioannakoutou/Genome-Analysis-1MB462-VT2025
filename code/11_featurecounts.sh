#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 10:00:00
#SBATCH -J featurecounts
#SBATCH --mail-type=ALL
#SBATCH --mail-user ioanna-kyriaki.koutou.9176@student.uu.se

#modules
module load bioinfo-tools
module load subread/2.0.3

#
export GFF_FILE="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/02_genome_annotation/eggnogmapper_output/eggnog_anot.emapper.decorated.gff"
export BAM_DIR="/proj/uppmax2025-3-3/nobackup/work/ioko9176/hisat2_output"
export OUTPUT_DIR="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/03_dea/featurecounts_output"

#run featureCounts
featureCounts -T 4 -p --countReadPairs \
  -t gene -g ID \
  -a "$GFF_FILE" \
  -o "${OUTPUT_DIR}/counts.txt" \
  -f \
  ${BAM_DIR}/*.bam
