#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=03:00:00
#SBATCH --job-name=generate_count_table
#SBATCH --mail-user=nina.eldridge@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/neldridge/rna_seq/data/10_generate_count_table/generate_count_table_%j.o
#SBATCH --error=/data/users/neldridge/rna_seq/fails/10_generate_count_table_%j.e
#SBATCH --partition=pall

#load required models
module load vital-it/latest
module load UHTS/Analysis/bedops/2.4.40;
module add UHTS/Analysis/subread/2.0.1

#define variable for main directory
main_dir=/data/users/neldridge/rna_seq
datadir=${main_dir}/data/10_generate_count_table
homo_dir=${main_dir}/data/3_prepare_annotation

cd $datadir

# counts the reads on CDS
featureCounts -T 4 -t CDS -g gene_id \
-a ${homo_dir}/GRCh38.p13.genome.gtf.gz \
-o ${datadir}/CDS_counts_rawfile.txt ${main_dir}/data/7_mapping_genome/*_GRCh38.p13_sorted.bam

# extract those reads that are mapped to different biotypes
featureCounts -T 4 -t exon -g gene_biotype \
-a ${homo_dir}/GRCh38.p13.genome.gtf.gz \
-o ${datadir}/biotype_counts_rawfile.txt ${main_dir}/data/7_mapping_genome/*_GRCh38.p13_sorted.bam


# extracts the gene ID and sample columns
cut -f 1,7-10 CDS_counts_rawfile.txt > CDS_counts_processed.txt

# extract biotypes and sample columns
cut -f 1,7-10 biotype_counts_rawfile.txt > biotype_counts_processed.txt