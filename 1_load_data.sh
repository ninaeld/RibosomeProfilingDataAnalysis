#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=get_data
#SBATCH --mail-user=nina.eldridge@students.unibe.ch
#SBATCH --mail-type=begin,fail,end
#SBATCH --output=/data/users/neldridge/rna_seq/data/1_raw_data/1_load_data_%j.o
#SBATCH --error=/data/users/neldridge/rna_seq/fails/1_load_data_%j.e
#SBATCH --partition=pall

#define variable for main directory
main_dir=/data/users/neldridge/rna_seq

#load modules
module load vital-it/latest
module load UHTS/Analysis/sratoolkit/2.10.7;

#vdb-config --interactive

#change directory to where data will be stored
cd $main_dir/data/1_raw_data

#download samples with SRA prefetch
prefetch --option-file samples.txt

#change sra files to zipped fastq files
fastq-dump --gzip SRR*.sra

