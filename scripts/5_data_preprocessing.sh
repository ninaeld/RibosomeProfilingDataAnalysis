#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=02:00:00
#SBATCH --job-name=data_preprocessing
#SBATCH --mail-user=nina.eldridge@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/neldridge/rna_seq/data/5_pre_processing/pre_processing_%j.o
#SBATCH --error=/data/users/neldridge/rna_seq/fails/5_pre_processing_%j.e
#SBATCH --partition=pall

#load required models
module load vital-it/latest
module load UHTS/Quality_control/cutadapt/2.5
module load UHTS/Quality_control/fastqc/0.11.9
module load UHTS/Analysis/MultiQC/1.8

#define variable for main directory
main_dir=/data/users/neldridge/rna_seq
datadir=${main_dir}/data/1_raw_data/human_fastq
outdir=${main_dir}/data/5_pre_processing

#singularity run-help /software/singularity/containers/Cutadapt-3.4-1.ubuntu20.sif;

cd $datadir

#iterate over all fastq files and cut adapters
for filename in *.gz; do
cutadapt ${filename} -j 4 -a CTGTAGGCACCATCAAT \
--discard-untrimmed \
-q 25 --minimum-length 25 \
-o ${filename%.fastq.gz}_clipped.fastq.gz;
done

# cut 4 nucleotides
for filename in *_clipped.fastq.gz; do
cutadapt ${filename} -j 4 -q 25 \
--cut -4 \
--minimum-length 25 \
-o ${filename%.fastq.gz}_trimmed.fastq.gz;
done

# do fastqc on clipped and trimmed reads
fastqc -o ${main_dir}/data/5_pre_processing -t 6 ${main_dir}/data/1_raw_data/human_fastq/*_trimmed.fastq.gz

# clip 3' adapter
singularity exec \
--bind ${main_dir} \
/software/singularity/containers/Cutadapt-3.4-1.ubuntu20.sif \
for x in $(ls -d *.fastq.gz); do echo ${x}; \
cutadapt -j 4 -a CTGTAGGCACCATCAAT \
-q 25 --minimum-length 25 --discard-untrimmed \
-o $(basename ${x} .fastq.gz)_clipped.fastq.gz \
${x} 1> $(basename ${x} .fastq.gz)_clipped_log.txt; done



