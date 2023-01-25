# 10) Differential expression analysis
# use DESeq2
# Nina Eldridge 23.12.22

setwd("C:/Users/Nina Eldridge/Documents/Master Bioinformatik/HS 2022/RNAseq/10_counts_table")

library(DESeq2)
#BiocManager::install("DESeq2")
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(BiocParallel)
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#install.packages("org.Rn.eg.db")
BiocManager::install("org.Hs.eg.db")
#BiocManager::install("AnnotationDbi")
library(org.Hs.eg.db) # Rat's Annotationdbi
register(MulticoreParam(4)) # Change this based on your computer core count

# Setting up color profiles from colorbrewer
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mycols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

sample_1 <- "WT"
sample_2 <- "KO"
no_of_reps <- 2

project_name <- paste("RPF", sample_2, "vs", sample_1, sep = "_")


## Create Samples dataframe
condition_column <- c(rep(sample_1, no_of_reps), rep(sample_2, no_of_reps))

run_column <- c(paste(sample_1, "1", sep = "_"), paste(sample_1, "2", sep = "_"),
                paste(sample_2, "1", sep = "_"), paste(sample_2, "2", sep = "_"))

rep_column <- c("A", "B", "A", "B")

samples_df <- data.frame(condition_column, run_column, rep_column)

colnames(samples_df) <- c("condition", "run", "rep")

rownames(samples_df) <- samples_df$run

samples_df$condition <- factor(rep(c(rep(sample_1, no_of_reps),
                                     rep(sample_2, no_of_reps))))


featurecount_data <- read.table("CDS_counts_processed.txt", header = TRUE, row.names = 1)


# Reorder to make the order consistent with samples$run
featurecount_data <- featurecount_data[ , c(3, 4, 1, 2)]

# Change colnames
colnames(featurecount_data) <- rownames(samples_df)

# Import as DESeqDataSet (dds)
dds <- DESeqDataSetFromMatrix(countData = featurecount_data, colData = samples_df, design = ~ condition)

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Factor levels
# WT columns should be first followed by KO/treatment
dds$condition <- factor(c(rep(sample_1, no_of_reps),
                          rep(sample_2, no_of_reps)),
                        levels = c(sample_1,
                                   sample_2))

# Differential expression analysis
dds <- DESeq(dds)



#### QC ####

# Log transformation for data quality assessment
rld <- rlog(dds, blind = FALSE)

# Sample distance matrix
sampleDists <- as.matrix(dist(t(assay(rld))))
pdf(paste(project_name, "QC_sample_distance_matrix_CDS.pdf", sep = "_"))
heatmap.2(as.matrix(sampleDists),
          key = TRUE,
          trace = "none",
          col = colorpanel(100, "#2b8cbe", "#7bccc4", "white"),
          ColSideColors = mycols[dds$condition],
          RowSideColors = mycols[dds$condition],
          margin = c(10, 10), main = "Sample Distance Matrix")
dev.off()


# Count matrix heatmap
select <- order(rowMeans(counts(dds, normalized = TRUE)))
df <- as.data.frame(colData(dds)[ , c("condition","rep")])

pdf(paste(project_name, "QC_count_matrix_CDS.pdf", sep = "_"))
pheatmap(assay(rld)[select, ],
         cluster_rows = FALSE,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         annotation_col = df)
dev.off()


# PCA plot

pdf(paste(project_name, "QC_PCA_CDS.pdf", sep = "_"))
pcaData <- plotPCA(rld, intgroup = c("condition", "rep"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData,
      aes(PC1, PC2, color = condition, shape = rep)) +
      geom_point(size = 3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      coord_fixed()

dev.off()


#### Resume the analysis ####

res_2_vs_1 <- results(dds, contrast = c("condition", sample_2, sample_1), alpha = 0.05)


convertIDs <- function( ids, from, to, db, ifMultiple = c("putNA", "useFirst")) {
  stopifnot(inherits(db, "AnnotationDb"))
  ifMultiple <- match.arg(ifMultiple)
  suppressWarnings(selRes <- AnnotationDbi::select(
    db, keys = ids, keytype = from, columns = c(from, to)))
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[duplicated(selRes[ , 1] ), 1]
    selRes <- selRes[!selRes[ , 1] %in% duplicatedIds, ]
  }
  return(selRes[match(ids, selRes[ , 1]), 2])
}

# # Check columns in the database:
columns(org.Hs.eg.db)

# Actual adding of the column

res_2_vs_1$GeneID <- row.names(res_2_vs_1)
res_2_vs_1$gene_symbol <- convertIDs(row.names(res_2_vs_1), "ENSEMBL", "SYMBOL", org.Hs.eg.db)


## Data summary
summary(res_2_vs_1)

## Convert DESeq object to dataframe
res_2_vs_1_df <- as.data.frame(res_2_vs_1)

## Define differentially expressed genes as:
## Upregulated (LFC > 0.5 & padj < 0.05)
## Downregulated (LFC < - 0.5 & padj < 0.05)

res_2_vs_1_df$regulation_level <- ifelse((res_2_vs_1_df$log2FoldChange > 0.5 & res_2_vs_1_df$padj < 0.05), "Upregulated",
                                         ifelse((res_2_vs_1_df$log2FoldChange < - 0.5 & res_2_vs_1_df$padj < 0.05), "Downregulated",
                                                "Unchanged"))

write.table (res_2_vs_1_df,
             file = paste(project_name, "DESeq2_res.csv", sep = "_"),
             sep = ",",
             row.names = FALSE,
             col.names = TRUE,
             quote = FALSE)



res_2_vs_1_df$regulation_level <- factor(res_2_vs_1_df$regulation_level, levels = c("Upregulated", "Downregulated", "Unchanged"))

## Remove rows with NA in padj column before plotting
res_2_vs_1_df <- res_2_vs_1_df[!is.na(res_2_vs_1_df$padj), ]

pdf(paste(project_name, "Volcano_plot.pdf", sep = "_"), width = 4, height = 5)
ggplot(res_2_vs_1_df,
      aes(x = -log10(padj), y = log2FoldChange, color = regulation_level)) +
      geom_point(size = 1) +
      scale_color_manual(values = c("#404788FF", "#73D055FF", "#999999")) +
      xlab("-log10(adjusted p-value)") +
      ylab("Log2 fold change") +
      labs(color = "Regulation level") +
      theme_bw()
dev.off()

