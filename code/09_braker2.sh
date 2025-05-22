#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 2-00:00:00
#SBATCH -J braker2_annot
#SBATCH --mail-type=ALL
#SBATCH --mail-user ioanna-kyriaki.koutou.9176@student.uu.se

#modules
module load bioinfo-tools
module load samtools/1.19
module load braker/2.1.6

#in_out paths
export GENOME="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/02_genome_annotation/repeatmasker_output/polished_assembly.fasta.masked"
export BAM_DIR="/proj/uppmax2025-3-3/nobackup/work/ioko9176/hisat2_output"
export BAM_IN="$BAM_DIR/merged_all_sorted.bam"
export PROTEIN="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/data/protein/embryophyte_proteomes.faa"
export BRAKER_OUT="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses/02_genome_annotation/braker2_output"

#environment paths
export AUGUSTUS_CONFIG_PATH="/home/ioko9176/augustus_config"
export AUGUSTUS_BIN_PATH="/sw/bioinfo/augustus/3.4.0/snowy/bin"
export AUGUSTUS_SCRIPTS_PATH="/sw/bioinfo/augustus/3.4.0/snowy/scripts"
export GENEMARK_PATH="/sw/bioinfo/GeneMark/4.68-es/snowy"

#genemark license
cp -vf /sw/bioinfo/GeneMark/keyfile/gm_key $HOME/.gm_key
ls -l "$HOME"/.gm_key


#step 1: merge bam files
cd "$BAM_DIR"
samtools merge -@ 8 "$BAM_IN" Control_*_sorted.bam Heat_treated_*_sorted.bam
samtools index "$BAM_IN"

#step 2: run braker
braker.pl \
    --genome="$GENOME" \
    --bam="$BAM_IN" \
    --prot_seq="$PROTEIN" \
    --etpmode \
    --softmasking \
    --species=niphotrichum_japonicum \
    --cores=8 \
    --gff3 \
    --workingdir="$BRAKER_OUT"
