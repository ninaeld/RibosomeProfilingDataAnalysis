#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=03:00:00
#SBATCH --job-name=remove_undesired
#SBATCH --mail-user=nina.eldridge@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/neldridge/rna_seq/data/6_remove_undesired_RNA/remove_RNA_%j.o
#SBATCH --error=/data/users/neldridge/rna_seq/fails/6_remove_undesired_RNA_%j.e
#SBATCH --partition=pall

#load required models
module load vital-it/latest
module load UHTS/Aligner/bowtie/1.2.0
module add UHTS/Analysis/samtools/1.10

#define variable for main directory
main_dir=/data/users/neldridge/rna_seq
datadir=${main_dir}/data/1_raw_data/human_fastq
RNAdir=${main_dir}/data/4_bowtie

cd $datadir

# map the reads to undesired RNA
# only keep those that are not mapped!
# do it with clipped and trimmed reads

# unzip reads
#for filename in *_trimmed.fastq.gz;do
#gunzip ${filename}
#done

# with bowtie align reads
for filename in *_trimmed.fastq;do
bowtie -S -t -p 4 \
$RNAdir/GRCh38_p13_r-t-sno-sn-RNA_ENSEMBL_NCBI_GtRNAdb  ${filename} \
--un "$main_dir/data/6_remove_undesired_RNA/${filename%.fastq}_desired_reads.fastq" \
2> "${filename%.fastq}_desired_reads_log.txt" > /dev/null;
done