#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -J soft_masking
#SBATCH --mail-type=ALL
#SBATCH --mail-user ioanna-kyriaki.koutou.9176@student.uu.se

#modules
module load bioinfo-tools
module load RepeatModeler/2.0.4
module load RepeatMasker/4.1.5

#
export WORKDIR="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/02_genome_annotation"
export IN="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/01_genome_assembly/pilon_output/polished_assembly.fasta"
export RMODELER_OUT="$WORKDIR/repeatmasker_output/repeatmodeler"
export RMASKER_OUT="$WORKDIR/repeatmasker_output"
export DB_NAME="njap_chr3_db"
export THREADS=4

#step 1: build the db
cd "$RMODELER_OUT"
BuildDatabase -name "$DB_NAME" -engine ncbi "$IN"

#step 2: run RepeatModeler (custom repeat library)
RepeatModeler -database "$DB_NAME" --threads "$THREADS"

#step 3: run RepeatMasker
cd $RMASKER_OUT
LIB=$(find "$RMODELER_OUT" -type f -name "consensi.fa.classified" | head -n 1)
RepeatMasker -pa "$THREADS" -lib "$LIB" -xsmall -gff -html -dir "$RMASKER_OUT" "$IN"
