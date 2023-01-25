setwd("/Users/puneetsharma/Documents/Work/Teaching/RNA_seq_2022/BAM_files/")

# Load the package
library(RiboseQC)

## Prepare Annotation

prepare_annotation_files(annotation_directory = "~/Documents/Work/Annotation/Humans/Genome/Test/",
                         twobit_file = "~/Documents/Work/Annotation/Humans/Genome/GRCh38.p13.genome.2bit",
                         gtf_file = "~/Documents/Work/Annotation/Humans/GTF/Homo_Sapiens_gencode_v35_annotation.gtf",
                         scientific_name = "Homo.sapiens",
                         annotation_name = "GRCh38_p13",
                         export_bed_tables_TxDb = F,
                         forge_BSgenome = T,
                         create_TxDb = T)


genome_bam <- c("/Users/puneetsharma/Documents/Work/Teaching/RNA_seq_2022/BAM_files/RPF_WT_Rep1_clpd_tr_no_r-t-sno-sn-RNA_GRCh38_p13_sorted.bam",
                "/Users/puneetsharma/Documents/Work/Teaching/RNA_seq_2022/BAM_files/RPF_WT_Rep2_clpd_tr_no_r-t-sno-sn-RNA_GRCh38_p13_sorted.bam",
                "/Users/puneetsharma/Documents/Work/Teaching/RNA_seq_2022/BAM_files/RPF_KO_Rep1_clpd_tr_no_r-t-sno-sn-RNA_GRCh38_p13_sorted.bam",
                "/Users/puneetsharma/Documents/Work/Teaching/RNA_seq_2022/BAM_files/RPF_KO_Rep2_clpd_tr_no_r-t-sno-sn-RNA_GRCh38_p13_sorted.bam")


load_annotation("~/Documents/Work/Annotation/Humans/Genome/Homo_Sapiens_gencode_v35_annotation.gtf_Rannot")

###### QC plots ######

RiboseQC_analysis(annotation_file ="~/Documents/Work/Annotation/Humans/Genome/Homo_Sapiens_gencode_v35_annotation.gtf_Rannot",
                  bam_files = genome_bam,
                  fast_mode = T,
                  report_file = "RPF_samples_QC.html",
                  sample_names = c("WT_Rep1", "WT_Rep2",
                                   "KO_Rep1", "KO_Rep2"),
                  dest_names = c("WT_Rep1", "WT_Rep2",
                                 "KO_Rep1", "KO_Rep2"),
                  write_tmp_files = F)
