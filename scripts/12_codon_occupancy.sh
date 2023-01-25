#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=03:00:00
#SBATCH --job-name=codon_occupancy
#SBATCH --mail-user=nina.eldridge@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/neldridge/rna_seq/data/12_codon_occupancy/codon_occupancy_%j.o
#SBATCH --error=/data/users/neldridge/rna_seq/fails/12_codon_occupancy_%j.e
#SBATCH --partition=pall

#load the required module
module load UHTS/Aligner/bowtie/1.2.0

#define variable for main directory
main_dir=/data/users/neldridge/rna_seq

#fast q files in 6 desired reads
cd $main_dir/data/6_remove_undesired_RNA


for filename in *_desired_reads.fastq;do
bowtie -t -p 4 -v 1 -m 1 \
--best \
--strata \
--norc \
$main_dir/data/4_bowtie/GRCh38_p13_APPRIS_CDS_plus18 \
-q ${filename} \
-S "$main_dir/data/12_codon_occupancy/${filename%.fastq}_GRCh38_p13_APPRIS_CDS.sam" \
2> "$main_dir/data/12_codon_occupancy/${filename%.fastq}_GRCh38_p13_APPRIS_CDS_log.txt";
done
