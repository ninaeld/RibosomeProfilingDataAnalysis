#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=02:00:00
#SBATCH --job-name=multiqc
#SBATCH --mail-user=nina.eldridge@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/neldridge/rna_seq/data/multiqc_%j.o
#SBATCH --error=/data/users/neldridge/rna_seq/fails/multiqc_%j.e
#SBATCH --partition=pall

#define variable for main directory
main_dir=/data/users/neldridge/rna_seq

# first the fastQC module has to be loaded
module load UHTS/Analysis/MultiQC/1.8

multiqc .