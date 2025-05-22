#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 09:00:00
#SBATCH -J hisat2_alignment
#SBATCH --mail-type=ALL
#SBATCH --mail-user ioanna-kyriaki.koutou.9176@student.uu.se

#modules
module load bioinfo-tools
module load samtools/1.20
module load HISAT2/2.2.1

#
export READ_DIR="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/data/rna_seq"
export GENOME_FA="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/02_genome_annotation/repeatmasker_output/polished_assembly.fasta.masked"
export OUT="/proj/uppmax2025-3-3/nobackup/work/ioko9176/hisat2_output"
export GENOME_INDEX="/proj/uppmax2025-3-3/nobackup/work/ioko9176/genome_index/masked_polished_index"
export THREADS="4"

#step 1: index the masked genome
hisat2-build "$GENOME_FA" "$GENOME_INDEX"

#step 2: run hisat2 
for FWD in "$READ_DIR"/*_f1.fq.gz; do
    SAMPLE=$(basename "$FWD" _f1.fq.gz)
    
    #reverse read filename
    REV="$READ_DIR/${SAMPLE}_r2.fq.gz"

    hisat2 -x "$GENOME_INDEX" -1 "$FWD" -2 "$REV" --threads "$THREADS" \
    | samtools view -bS -\
    | samtools sort -@ "$THREADS" -o "$OUT/${SAMPLE}_sorted.bam"

    samtools index "$OUT/${SAMPLE}_sorted.bam"
done
 
# symbolic links back to the home dir             
ln -s "$OUT"/*.bam /home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/02_genome_annotation/hisat2_output/
ln -s "$OUT"/*.bai /home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/02_genome_annotation/hisat2_output/
ln -s "$GENOME_INDEX".* /home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/02_genome_annotation/hisat2_output/genome_index/


