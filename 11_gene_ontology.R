# 11) Gene ontology
# Nina Eldridge 4.1.2023

setwd("C:/Users/Nina Eldridge/Documents/Master Bioinformatik/HS 2022/RNAseq/11_gene_ontology")

library(grid)
library(gridExtra)
#BiocManager::install("pathview")
library(pathview)
library(AnnotationDbi)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
#BiocManager::install("topGO")
library(topGO)
library(ggplot2)

sample_name = "RPF_KO_vs_WT"

## Load DESeq2 output file.

df <- read.csv("RPF_KO_vs_WT_DESeq2_res.csv", sep = ",", header = TRUE)

## Set rownames as GeneID column

rownames(df) <- df$GeneID
df <- df[order(df$padj), ]

## Define Up and Down regulated genes
## Here we are considering padj values only so that we do not lose important terms

genes_up <- which(df$padj < 0.05 & df$log2FoldChange > 0)
genes_down <- which(df$padj < 0.05 & df$log2FoldChange < 0)

all_genes_names <- rownames(df)

genes_up <- rownames(df)[genes_up]
genes_down <- rownames(df)[genes_down]

genelist_up <- factor(as.integer(all_genes_names %in% genes_up))
names(genelist_up) <- all_genes_names

genelist_down <- factor(as.integer(all_genes_names %in% genes_down))
names(genelist_down) <- all_genes_names

## Make sure to change your "mapping" DB to your organism of interest
allGO2genes <- annFUN.org(whichOnto = "ALL",
                          feasibleGenes = NULL,
                          mapping = "org.Hs.eg.db",
                          ID = "ensembl")

# Create TopGOData object (Biological processes; Upregulated)
GOdata_up_bp <- new("topGOdata",
                    ontology = "BP",
                    allGenes = genelist_up,
                    annot = annFUN.GO2genes,
                    GO2genes = allGO2genes,
                    nodeSize = 10)
# Create TopGOData object (Molecular function; Upregulated)
GOdata_up_mf <- new("topGOdata",
                    ontology = "MF",
                    allGenes = genelist_up,
                    annot = annFUN.GO2genes,
                    GO2genes = allGO2genes,
                    nodeSize = 10)

# Create TopGOData object (Cellular components; Upregulated)
GOdata_up_cc <- new("topGOdata",
                    ontology = "CC",
                    allGenes = genelist_up,
                    annot = annFUN.GO2genes, 
                    GO2genes = allGO2genes,
                    nodeSize = 10)
# Create TopGOData object (Biological processes; Downregulated)
GOdata_down_bp <- new("topGOdata",
                      ontology = "BP",
                      allGenes = genelist_down,
                      annot = annFUN.GO2genes,
                      GO2genes = allGO2genes,
                      nodeSize = 10)

# Create TopGOData object (Molecular function; Downregulated)
GOdata_down_mf <- new("topGOdata",
                      ontology = "MF",
                      allGenes = genelist_down,
                      annot = annFUN.GO2genes,
                      GO2genes = allGO2genes,
                      nodeSize = 10)

# Create TopGOData object (Cellular components; Downregulated)
GOdata_down_cc <- new("topGOdata",
                      ontology = "CC",
                      allGenes = genelist_down,
                      annot = annFUN.GO2genes,
                      GO2genes = allGO2genes,
                      nodeSize = 10)

resultFis_up_bp <- runTest(GOdata_up_bp, algorithm = "elim", statistic = "fisher")
resultFis_up_mf <- runTest(GOdata_up_mf, algorithm = "elim", statistic = "fisher")
resultFis_up_cc <- runTest(GOdata_up_cc, algorithm = "elim", statistic = "fisher")
resultFis_down_bp <- runTest(GOdata_down_bp, algorithm = "elim", statistic = "fisher")
resultFis_down_mf <- runTest(GOdata_down_mf, algorithm = "elim", statistic = "fisher")
resultFis_down_cc <- runTest(GOdata_down_cc, algorithm = "elim", statistic = "fisher")


# Function to extract data from S4 objects

parse_tables <- function(GO_data, statistics)
{
  goEnrichment <- GenTable(GO_data, weightFisher = statistics, topNodes = 20)
  sub("< ", "", goEnrichment$weightFisher)
  goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep = ", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels = rev(goEnrichment$Term))
  goEnrichment$weightFisher <- as.numeric(sub("< ", "", goEnrichment$weightFisher))  
  goEnrichment
}

GOres_up_bp <- parse_tables(GOdata_up_bp, resultFis_up_bp)
GOres_up_mf <- parse_tables(GOdata_up_mf, resultFis_up_mf)
GOres_up_cc <- parse_tables(GOdata_up_cc, resultFis_up_cc)

GOres_down_bp <- parse_tables(GOdata_down_bp, resultFis_down_bp)
GOres_down_mf <- parse_tables(GOdata_down_mf, resultFis_down_mf)
GOres_down_cc <- parse_tables(GOdata_down_cc, resultFis_down_cc)


# Funciton for the plots

plot_GO <- function(GO_data, Ontology, Regulation, use_color) {
  GO_data$log_weightFisher <- (- log10(as.numeric(GO_data$weightFisher)))
  ggplot(GO_data, 
         aes(x = GO_data$log_weightFisher,
             y = GO_data$Term)) +
    geom_segment(aes(x = 0,
                     xend = GO_data$log_weightFisher,
                     y = GO_data$Term,
                     yend = GO_data$Term),
                 colour = use_color)  +
    geom_point(aes(size = GO_data$Significant),
               colour = use_color) +
    scale_size_area(name = "Gene counts") +
    xlab("Enrichment (- log10 Pvalue)") +
    ylab(Ontology) +
    ggtitle(Regulation) +
    scale_x_continuous() +
    theme_bw() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"))
}

# Make Plots

plot_up_BP <- plot_GO(GOres_up_bp, "Biological Proccess", "TopGO Up (fisher's exact test)", "#404788FF")
plot_up_MF <- plot_GO(GOres_up_mf, "Molecular Function", "TopGO Up (fisher's exact test)", "#404788FF")
plot_up_CC <- plot_GO(GOres_up_cc, "Cellular Component", "TopGO Up (fisher's exact test)", "#404788FF")


plot_down_BP <- plot_GO(GOres_down_bp, "Biological Proccess", "TopGO Down (fisher's exact test)", "#73D055FF")
plot_down_MF <- plot_GO(GOres_down_mf, "Molecular Function", "TopGO Down (fisher's exact test)", "#73D055FF")
plot_down_CC <- plot_GO(GOres_down_cc, "Cellular Component", "TopGO Down (fisher's exact test)", "#73D055FF")


# Export as PDF

pdf(paste(sample_name, "Biological_Proccess_TopGO_Up_fisher.pdf", sep = "_"))
plot_up_BP
dev.off()

pdf(paste(sample_name, "Molecular_Function_TopGO_Up_fisher.pdf", sep = "_"))
plot_up_MF
dev.off()

pdf(paste(sample_name, "Cellular_Component_TopGO_Up_fisher.pdf", sep = "_"))
plot_up_CC
dev.off()

pdf(paste(sample_name, "Biological_Proccess_TopGO_down_fisher.pdf", sep = "_"))
plot_down_BP
dev.off()

pdf(paste(sample_name, "Molecular_Function_TopGO_down_fisher.pdf", sep = "_"))
plot_down_MF
dev.off()

pdf(paste(sample_name, "Cellular_Component_TopGO_down_fisher.pdf", sep = "_"))
plot_down_CC
dev.off()



# retrieve genes2GO list from the "expanded" annotation in GOdata
# SOURCE: https://support.bioconductor.org/p/29775/

# For upregulated
all_GO_up_bp = genesInTerm(GOdata_up_bp)
all_GO_up_mf = genesInTerm(GOdata_up_mf)
all_GO_up_cc = genesInTerm(GOdata_up_cc)
# For down regulated
all_GO_down_bp = genesInTerm(GOdata_down_bp)
all_GO_down_mf = genesInTerm(GOdata_down_mf)
all_GO_down_cc = genesInTerm(GOdata_down_cc)

data_to_extract_up_bp = lapply(all_GO_up_bp, function(x) x[x %in% genes_up])
data_to_extract_up_mf = lapply(all_GO_up_mf, function(x) x[x %in% genes_up])
data_to_extract_up_cc = lapply(all_GO_up_cc, function(x) x[x %in% genes_up])

data_to_extract_down_bp = lapply(all_GO_down_bp, function(x) x[x %in% genes_down])
data_to_extract_down_mf = lapply(all_GO_down_mf, function(x) x[x %in% genes_down])
data_to_extract_down_cc = lapply(all_GO_down_cc, function(x) x[x %in% genes_down])



# To view your significant genes for, say, GO:0051427
# data_to_extract_up_bp[["GO:0032543"]]

# Extract our genes in all identified GO terms:
# Source: https://stackoverflow.com/a/15201690

no_of_observations_up_bp <- sapply(data_to_extract_up_bp, length)
no_of_observations_up_mf <- sapply(data_to_extract_up_mf, length)
no_of_observations_up_cc <- sapply(data_to_extract_up_cc, length)

no_of_observations_down_bp <- sapply(data_to_extract_down_bp, length)
no_of_observations_down_mf <- sapply(data_to_extract_down_mf, length)
no_of_observations_down_cc <- sapply(data_to_extract_down_cc, length)

seq_max_up_bp <- seq_len(max(no_of_observations_up_bp))
seq_max_up_mf <- seq_len(max(no_of_observations_up_mf))
seq_max_up_cc <- seq_len(max(no_of_observations_up_cc))

seq_max_down_bp <- seq_len(max(no_of_observations_down_bp))
seq_max_down_mf <- seq_len(max(no_of_observations_down_mf))
seq_max_down_cc <- seq_len(max(no_of_observations_down_cc))


GO_df_up_bp <- t(sapply(data_to_extract_up_bp, "[", i = seq_max_up_bp))
GO_df_up_mf <- t(sapply(data_to_extract_up_mf, "[", i = seq_max_up_mf))
GO_df_up_cc <- t(sapply(data_to_extract_up_cc, "[", i = seq_max_up_cc))

GO_df_down_bp <- t(sapply(data_to_extract_down_bp, "[", i = seq_max_down_bp))
GO_df_down_mf <- t(sapply(data_to_extract_down_mf, "[", i = seq_max_down_mf))
GO_df_down_cc <- t(sapply(data_to_extract_down_cc, "[", i = seq_max_down_cc))


write.table(GO_df_up_bp,
            file = paste(sample_name, "BP_Up_fisher_Genes_in_GO.txt", sep = "_"),
            sep = "\t",
            row.names = TRUE,
            col.names = FALSE,
            quote = FALSE,
            na = "")

write.table(GO_df_up_mf,
            file = paste(sample_name, "MF_Up_fisher_Genes_in_GO.txt", sep = "_"),
            sep = "\t",
            row.names = TRUE,
            col.names = FALSE,
            quote = FALSE,
            na = "")

write.table(GO_df_up_cc,
            file = paste(sample_name, "CC_Up_fisher_Genes_in_GO.txt", sep = "_"),
            sep = "\t",
            row.names = TRUE,
            col.names = FALSE,
            quote = FALSE,
            na = "")

write.table(GO_df_down_bp,
            file = paste(sample_name, "BP_Down_fisher_Genes_in_GO.txt", sep = "_"),
            sep = "\t",
            row.names = TRUE,
            col.names = FALSE,
            quote = FALSE,
            na = "")

write.table(GO_df_down_mf,
            file = paste(sample_name, "MF_Down_fisher_Genes_in_GO.txt", sep = "_"),
            sep = "\t",
            row.names = TRUE,
            col.names = FALSE,
            quote = FALSE,
            na = "")

write.table(GO_df_down_cc,
            file = paste(sample_name, "CC_Down_fisher_Genes_in_GO.txt", sep = "_"),
            sep = "\t",
            row.names = TRUE,
            col.names = FALSE,
            quote = FALSE,
            na = "")
