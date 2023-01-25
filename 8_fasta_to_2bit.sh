#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=01:00:00
#SBATCH --job-name=fasta_to_2bit
#SBATCH --mail-user=nina.eldridge@students.unibe.ch
#SBATCH --mail-type=begin,fail,end
#SBATCH --output=/data/users/neldridge/rna_seq/data/8_fasta_to_2bit/to_2bit_%j.o
#SBATCH --error=/data/users/neldridge/rna_seq/fails/8_to_2bit_%j.e
#SBATCH --partition=pall

#load required models
module load vital-it/latest
module load SequenceAnalysis/blat/36 

main_dir=/data/users/neldridge/rna_seq
datadir=${main_dir}/data/3_prepare_annotation

cd $datadir

faToTwoBit GRCh38.p13.genome.fa ${main_dir}/data/8_fasta_to_2bit/GRCh38.p13.genome.2bit