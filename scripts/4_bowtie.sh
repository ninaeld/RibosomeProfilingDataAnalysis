#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --job-name=4_bowtie
#SBATCH --mail-user=nina.eldridge@students.unibe.ch
#SBATCH --mail-type=begin,fail,end
#SBATCH --output=/data/users/neldridge/rna_seq/data/4_bowtie/bowtie_%j.o
#SBATCH --error=/data/users/neldridge/rna_seq/data/4_bowtie/bowtie_%j.e
#SBATCH --partition=pall

module load UHTS/Aligner/bowtie/1.2.0;

main_dir=/data/users/neldridge/rna_seq
workdir=$main_dir/data/4_bowtie
datadir=$main_dir/data/3_prepare_annotation

cd $workdir

# For the "undesired" RNAs
bowtie-build $datadir/GRCh38_p13_r-t-sno-sn-RNA_ENSEMBL_NCBI_GtRNAdb.fa GRCh38_p13_r-t-sno-sn-RNA_ENSEMBL_NCBI_GtRNAdb

# For the genome
bowtie-build $datadir/GRCh38.p13.genome.fa GRCh38.p13.genome

# For the transcriptome
bowtie-build $datadir/GRCh38_p13_APPRIS_CDS_plus18.fa GRCh38_p13_APPRIS_CDS_plus18


# process transriptome
cd $datadir
# transcriptome FASTA is in multiline format, change to single line format
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' \
< GRCh38_p13_APPRIS_CDS_plus18.fa \
> GRCh38_p13_APPRIS_CDS_plus18_SingleLine.fa