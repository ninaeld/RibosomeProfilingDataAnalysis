./Codon_occupancy_cal.sh \
/data/users/neldridge/rna_seq/data/3_prepare_annotation/GRCh38_p13_APPRIS_CDS_plus18_SingleLine.fa \
/data/users/neldridge/rna_seq/data/12_codon_occupancy/RPF_KO_Rep1_clipped_trimmed_desired_reads_GRCh38_p13_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_KO_Rep1_Codon_occupancy.txt

./Codon_occupancy_cal.sh \
/data/users/neldridge/rna_seq/data/3_prepare_annotation/GRCh38_p13_APPRIS_CDS_plus18_SingleLine.fa \
/data/users/neldridge/rna_seq/data/12_codon_occupancy/RPF_KO_Rep2_clipped_trimmed_desired_reads_GRCh38_p13_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_KO_Rep2_Codon_occupancy.txt

./Codon_occupancy_cal.sh \
/data/users/neldridge/rna_seq/data/3_prepare_annotation/GRCh38_p13_APPRIS_CDS_plus18_SingleLine.fa \
/data/users/neldridge/rna_seq/data/12_codon_occupancy/RPF_WT_Rep1_clipped_trimmed_desired_reads_GRCh38_p13_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_WT_Rep1_Codon_occupancy.txt

./Codon_occupancy_cal.sh \
/data/users/neldridge/rna_seq/data/3_prepare_annotation/GRCh38_p13_APPRIS_CDS_plus18_SingleLine.fa \
/data/users/neldridge/rna_seq/data/12_codon_occupancy/RPF_WT_Rep2_clipped_trimmed_desired_reads_GRCh38_p13_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_WT_Rep2_Codon_occupancy.txt