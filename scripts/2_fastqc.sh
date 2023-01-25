#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=02:00:00
#SBATCH --job-name=fastqc
#SBATCH --mail-user=nina.eldridge@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/neldridge/rna_seq/data/2_fastqc/2_fastqc_%j.o
#SBATCH --error=/data/users/neldridge/rna_seq/fails/2_fast_qc_%j.e
#SBATCH --partition=pall

#define variable for main directory
main_dir=/data/users/neldridge/rna_seq

# first the fastQC module has to be loaded
module add UHTS/Quality_control/fastqc/0.11.9;

# run fastqc for all six reads
fastqc -o ${main_dir}/data/2_fastqc -t 6 ${main_dir}/data/1_raw_data/human_fastq/*.fastq.gz