#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 1-00:00:00
#SBATCH -J Pilon_polish
#SBATCH --mail-type=ALL
#SBATCH --mail-user ioanna-kyriaki.koutou.9176@student.uu.se

#modules
module load bioinfo-tools
module load bwa/0.7.18
module load samtools/1.20
module load Pilon/1.24

#working dir
export WORKDIR="/home/ioko9176/genome_analysis/Genome-Analysis-1MB462-VT2025/analyses"

#read paths
export READ1="$WORKDIR/00_read_preprocessing/trimmomatic_PE_output/chr3_illumina_R1_trimmed_paired.fastq.gz"
export READ2="$WORKDIR/00_read_preprocessing/trimmomatic_PE_output/chr3_illumina_R2_trimmed_paired.fastq.gz"

#output dir
export OUTDIR="$WORKDIR/01_genome_assembly/pilon_output"

#dir for intermediate files
export TMPDIR="$OUTDIR/intermediate_files"
mkdir -p "$TMPDIR"

#create index dir and symlink FASTA
mkdir -p "$OUTDIR/index_output"
ln -s "$WORKDIR/01_genome_assembly/flye_output/assembly.fasta" "$OUTDIR/index_output/assembly.fasta"

#genome path
export GENOME="$OUTDIR/index_output/assembly.fasta"

#step 1: index the genome
bwa index "$GENOME"

#step 2 & 3: map reads and produce a sorted BAM file
bwa mem -t 4 "$GENOME" "$READ1" "$READ2" | samtools view -bS - | samtools sort -o "$TMPDIR/pilon_alignment.sorted.bam"

#index the sorted BAM file
samtools index "$TMPDIR/pilon_alignment.sorted.bam"

#step 4: run pilon - polish the genome
java -Xmx16G -jar "$PILON_HOME/pilon.jar" \
  --genome "$GENOME" \
  --frags "$TMPDIR/pilon_alignment.sorted.bam" \
  --output "$OUTDIR/polished_assembly" \
  --threads 4 

