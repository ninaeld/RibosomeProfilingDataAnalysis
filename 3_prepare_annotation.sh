#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --job-name=3_prepare_annotation
#SBATCH --mail-user=nina.eldridge@students.unibe.ch
#SBATCH --mail-type=begin,fail,end
#SBATCH --output=/data/users/neldridge/rna_seq/data/3_prepare_annotation/catting_%j.o
#SBATCH --error=/data/users/neldridge/rna_seq/data/3_prepare_annotation/catting_%j.e
#SBATCH --partition=pall

main_dir=/data/users/neldridge/rna_seq
workdir=$main_dir/data/3_prepare_annotation

#change to working directory
cd $workdir

#concatenate all text files into one fasta file
cat *.txt > GRCh38_p13_r-t-sno-sn-RNA_ENSEMBL_NCBI_GtRNAdb.fa