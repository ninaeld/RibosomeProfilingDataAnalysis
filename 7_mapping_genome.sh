#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=03:00:00
#SBATCH --job-name=mapping_genome
#SBATCH --mail-user=nina.eldridge@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/neldridge/rna_seq/data/7_mapping_genome/mapping_genome_%j.o
#SBATCH --error=/data/users/neldridge/rna_seq/fails/7_mapping_genome_%j.e
#SBATCH --partition=pall

#load required models
module load vital-it/latest
module load UHTS/Aligner/bowtie/1.2.0
module add UHTS/Analysis/samtools/1.10

#define variable for main directory
main_dir=/data/users/neldridge/rna_seq
datadir=${main_dir}/data/1_raw_data/human_fastq
desiredRNAdir=${main_dir}/data/6_remove_undesired_RNA
mapGenomedir=/${main_dir}/data/4_bowtie

cd $desiredRNAdir


# go through the 4 fastq files with the desired RNA and map to genome
for filename in *_desired_reads.fastq;do

#--best reports alignments from best to worst
# --strata report only the alignments in best alignment stratum ->
# having the least nb of mismatches -> if strata, then also best specified
# -v mode no more than v mismatches
# strata in -v alignment mode: alignments stratum defined as total nb of mismathces
# in entire alignment
# -m --best and --strata are reporting modes -> -m 1 refrain from reporting
# any alignments for reads having more than 1 reportable alignment
# having all three means principled but weaker form of uniqueness
# Output -t print amount of time taken by each phase
# SAM S print in SAM format
# Performance -p amount of threads
bowtie -S -t -p 4 -v 1 -m 1 --best --strata \
$mapGenomedir/GRCh38.p13.genome \
-q ${filename} \
2> ${filename%.fastq}_GRCh38.p13.genome_log.txt | \
samtools view -h -F 4 -b > ${main_dir}/data/7_mapping_genome/"${filename%.fastq}_GRCh38.p13.bam";
done

# change to where bam files are stored
cd $main_dir/data/7_mapping_genome

# sort each bam file
for filename in *.bam;do
samtools sort -@ 4 $filename -o ${filename%.bam}_sorted.bam;
done

# remove unsorted bam file
# rm *GRCh38_p13.bam