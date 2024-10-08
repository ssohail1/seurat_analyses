#### with Kallisto ####
# libraries
library(Seurat)
library(Matrix)
library(tidyverse)
library(patchwork)
library(ggplot2)
# theme_set(theme_bw())

# # Install kb (includes installing kallisto and bustools)
# system("pip3 install kb-python", intern=TRUE)
# 
# # download an index
# system("kb ref -d mouse -i index.idx -g t2g.txt -f1 transcriptome.fasta",intern=TRUE)
# 
# # Pseudoalignment and counting
# ## system("kb count -i index.idx -g t2g.txt -x 10xv3 -o output --filter bustools -t 2 pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz",intern=TRUE)
# system("kb count -i index.idx -g t2g.txt -f1 /Users/sidrasohail/Downloads/lungCan_data_seurat/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz -x 10xv3 -o output --filter bustools -t 17 /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507943_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507943_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507945_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507945_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507947_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507947_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507948_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507948_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507949_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507949_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507950_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507950_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507951_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507951_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507952_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507952_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507953_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507953_S1_R2_001.fastq.gz  /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507954_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507954_S1_R2_001.fastq.gz  /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507955_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507955_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507956_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507956_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507957_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507957_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507958_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507958_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507959_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507959_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507960_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507960_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507961_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507961_S1_R2_001.fastq.gz",intern=TRUE)
# /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507943_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507943_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507945_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507945_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507947_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507947_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507948_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507948_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507949_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507949_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507950_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507950_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507951_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507951_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507952_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507952_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507953_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507953_S1_R2_001.fastq.gz  /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507954_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507954_S1_R2_001.fastq.gz  /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507955_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507955_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507956_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507956_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507957_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507957_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507958_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507958_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507959_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507959_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507960_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507960_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507961_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507961_S1_R2_001.fastq.gz

## first build the index with command 'kallisto index -i "name_of_output_indexfile.idx' [.idx is optional] FASTA-file(s)
#~/Downloads/kallisto/kallisto index -i mouseind.idx ~/Downloads/lungCan_data_seurat/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

## then run kallisto on paired or single-end FASTQ files
# kallisto quant -i indexfilename -o output pairA_1.fastq pairA_2.fastq pairB_1.fastq pairB_2.fastq
# ~/Downloads/kallisto/kallisto quant -i mouseind.idx -o output /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507943_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507943_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507945_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507945_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507947_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507947_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507948_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507948_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507949_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507949_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507950_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507950_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507951_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507951_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507952_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507952_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507953_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507953_S1_R2_001.fastq.gz  /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507954_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507954_S1_R2_001.fastq.gz  /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507955_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507955_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507956_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507956_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507957_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507957_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507958_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507958_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507959_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507959_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507960_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507960_S1_R2_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507961_S1_R1_001.fastq.gz /Users/sidrasohail/Downloads/lungCan_data_seurat/SRR19507961_S1_R2_001.fastq.gz


# Basic QC
list.files(".", recursive = TRUE)
# Slightly modified from BUSpaRse, just to avoid installing a few dependencies not used here
read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  ge <- ".genes.txt"
  genes <- readLines(file(paste0(dir, "/", name, ge)))
  barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}
res_mat <- read_count_output("./output/counts_unfiltered", name = "cells_x_genes")
dim(res_mat)
tot_counts <- colSums(res_mat)
lib_sat <- tibble(nCount = tot_counts, nGene = colSums(res_mat > 0))

options(repr.plot.width=9, repr.plot.height=6)
ggplot(lib_sat, aes(nCount, nGene)) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_x_log10() + scale_y_log10() + annotation_logticks()

ggplot(lib_sat, aes(nCount, nGene)) +
  geom_bin2d(bins = 50) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_x_log10() + scale_y_log10() + annotation_logticks()

summary(tot_counts)

#' @rdname knee_plot
#' @param mat Gene count matrix, a dgCMatrix.
#' @return `get_knee_df` returns a tibble with two columns: \code{total} for 
#' total UMI counts for each barcode, and \code{rank} for rank of the total 
#' counts, with number 1 for the barcode with the most counts.
#' @export
#' @importFrom dplyr row_number desc arrange
#' @importFrom Matrix colSums
get_knee_df <- function(mat) {
  total <- rank <- NULL
  tibble(total = Matrix::colSums(mat),
         rank = row_number(desc(total))) %>%
    distinct() %>%
    dplyr::filter(total > 0) %>% 
    arrange(rank)
}

#' @rdname knee_plot
#' @param df The data frame from \code{\link{get_knee_df}}.
#' @param lower Minimum total UMI counts for barcode for it to be considered
#' when calculating the inflection point; this helps to avoid the noisy part of
#' the curve for barcodes with very few counts.
#' @return `get_inflection` returns a \code{numeric(1)} for the total UMI count 
#' at the inflection point.
#' @note Code in part adapted from \code{barcodeRanks} from \code{DropetUtils}.
#' @export
#' @importFrom dplyr transmute
#' 
get_inflection <- function(df, lower = 100) {
  log_total <- log_rank <- total <-  NULL
  df_fit <- df %>% 
    dplyr::filter(total > lower) %>% 
    transmute(log_total = log10(total),
              log_rank = log10(rank))
  d1n <- diff(df_fit$log_total)/diff(df_fit$log_rank)
  right.edge <- which.min(d1n)
  10^(df_fit$log_total[right.edge])
}

#' Plot the transposed knee plot and inflection point
#' 
#' Plot a transposed knee plot, showing the inflection point and
#' the number of remaining cells after inflection point filtering. It's
#' transposed since it's more generalizable to multi-modal data. Taken from the 
#' BUSpaRse package.
#' 
#' @param df The data frame from \code{\link{get_knee_df}}.
#' @param inflection Output of \code{\link{get_inflection}}.
#' @return `knee_plot` returns a \code{ggplot2} object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_path geom_vline geom_hline 
#' scale_x_log10 scale_y_log10 labs annotation_logticks geom_text
knee_plot <- function(df, inflection) {
  total <- rank_cutoff <- NULL
  annot <- tibble(inflection = inflection,
                  rank_cutoff = max(df$rank[df$total > inflection]))
  ggplot(df, aes(total, rank)) +
    geom_path() +
    geom_vline(aes(xintercept = inflection), data = annot, linetype = 2, 
               color = "gray40") +
    geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2, 
               color = "gray40") +
    geom_text(aes(inflection, rank_cutoff, 
                  label = paste(rank_cutoff, "'cells'")),
              data = annot, vjust = 1) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = "Rank", x = "Total UMIs") +
    annotation_logticks()
}

options(repr.plot.width=9, repr.plot.height=6)
knee_df <- get_knee_df(res_mat)
inflection <- get_inflection(knee_df)
knee_plot(knee_df, inflection)

# Analysis

tr2g <- read_tsv("t2g.txt", col_names = c("transcript", "gene", "gene_name"))
tr2g <- distinct(tr2g[, c("gene", "gene_name")])

plot_pct_genes <- function(mat, tr2g, top_n = 20, symbol = "ensembl") {
  pct_tx <- rowSums(mat)
  gs <- rownames(mat)[order(-pct_tx)]
  df <- as.data.frame(t(mat[gs[1:20],]))
  df <- df %>%
    mutate_all(function(x) x/colSums(mat)) %>%
    pivot_longer(everything(), names_to = "gene")
  if (symbol == "ensembl") {
    df <- left_join(df, tr2g, by = "gene")
  } else {
    df <- rename(df, gene_name = gene)
  }
  df %>%
    mutate(gene = fct_reorder(gene_name, value, .fun = median)) %>%
    ggplot(aes(gene, value)) +
    geom_boxplot() +
    labs(x = "", y = "Proportion of total counts") +
    coord_flip()
}

options(repr.plot.width=6, repr.plot.height=10)
plot_pct_genes(res_mat, tr2g)

res_mat <- res_mat[, tot_counts > inflection]
res_mat <- res_mat[Matrix::rowSums(res_mat) > 0,]
dim(res_mat)

# Convert from Ensembl gene ID to gene symbol
rownames(res_mat) <- tr2g$gene_name[match(rownames(res_mat), tr2g$gene)]

(pbmc <- CreateSeuratObject(counts = res_mat, project = "pbmc1k", min.cells = 3, min.features = 200))

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# Visualize QC metrics as a violin plot
options(repr.plot.width=12, repr.plot.height=6)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 20)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# pbmc <- NormalizeData(pbmc)

options(repr.plot.width=6, repr.plot.height=10)
plot_pct_genes(GetAssayData(pbmc, slot = "counts"), tr2g, symbol = "symbol")
options(repr.plot.width=9, repr.plot.height=6)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc, log = FALSE)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
# When using repel, set xnudge and ynudge to 0 for optimal results


# centering and scaling the data matrix
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
# The scaling does not affect PCA or clustering results. 
# However, Seurat heatmaps (produced as shown below with DoHeatmap) 
# require genes in the heatmap to be scaled so that 
# highly-expressed genes don’t dominate. To make sure we don’t leave any genes 
# out of the heatmap later, we are scaling all genes in this tutorial.


# Principal component analysis
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Visualize
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
# which genes are contributing most to the first 2 PCs?
options(repr.plot.width=6, repr.plot.height=8)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

options(repr.plot.width=7, repr.plot.height=6)
FeaturePlot(pbmc, reduction = "pca", feature = "CST3")

# Determining dimensionality
options(repr.plot.width=9, repr.plot.height=6)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10, k.param = 20)
pbmc <- FindClusters(pbmc, resolution = 0.6)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# Run UMAP and t-SNE
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
options(repr.plot.width=16, repr.plot.height=5)
FeaturePlot(pbmc, reduction = "umap", features = c("CST3", "NKG7", "PPBP"),
            ncol = 3)


# Finding differentially expressed features (cluster biomarkers)
# Scanpy style gene rank plot
plot_gene_rank <- function(markers, n) {
  df_plot <- markers %>%
    group_by(cluster) %>%
    top_n(25, avg_log2FC) %>%
    mutate(rank = factor(row_number(desc(avg_log2FC))))
  ggplot(df_plot, aes(rank, avg_log2FC)) +
    geom_text(aes(label = gene), angle = -90, hjust = 1) +
    facet_wrap(~ cluster) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
}

# find markers for every cluster compared to all remaining cells, report only the positive ones
# default is wilcoxon
pbmc.markers <- FindAllMarkers(pbmc, test.use = "wilcox", only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
head(pbmc.markers)


options(repr.plot.width=12, repr.plot.height=8)
plot_gene_rank(pbmc.markers, 25)

# student's t-test
pbmc.markers.t <- FindAllMarkers(pbmc, test.use = "t", only.pos = TRUE, 
                                 min.pct = 0.25, logfc.threshold = 0.25)
plot_gene_rank(pbmc.markers.t, 25)


# Also logistic regression to test how good each gene is for deciding whether a cell is in a cluster
pbmc.markers.lr <- FindAllMarkers(pbmc, test.use = "LR", only.pos = TRUE, 
                                  min.pct = 0.25, logfc.threshold = 0.25)
plot_gene_rank(pbmc.markers.lr, 25)


marker_genes <- sort(c('IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                       'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',  
                       'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP', 'CCR7',
                       'S100A4'))

options(repr.plot.width=16, repr.plot.height=20)
VlnPlot(pbmc, features = marker_genes, ncol = 4)

options(repr.plot.width=16, repr.plot.height=20)
FeaturePlot(pbmc, features = marker_genes, ncol = 4)

# Assigning cell type identity to clusters
options(repr.plot.width=6, repr.plot.height=7)
DotPlot(pbmc, assay = "RNA", features = marker_genes, scale.by = "size") +
  coord_flip()

options(repr.plot.width=9, repr.plot.height=6)
new.cluster.ids <- c("CD14+ Mono", "Memory CD4 T", "Naive CD4 T", "B1", "FCGR3A+ Mono", 
                     "NK", "CD8+ T", "B2", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6) + NoLegend()


## replicate with full single-cell dataset
               
# Before Seurat - downloading data from Single-cell RNA sequencing reveals myeloid and T cell co-stimulation mediated by IL-7 anti-cancer immunotherapy study
# bioproject number is PRJNA844355; GEO: GSE205307

# got accession list from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP378029&o=acc_s%3Aa for sra files
# ran the following to retrieve sra data with prefetch:

# for sample in $(cat SRR_Acc_List.txt):
#   do
## ~/sratoolkit.3.0.7-mac64/bin/prefetch ${sample}
# prefetch ${sample}
# done

# for sample in $(cat SRR_Acc_List.txt):
# do
# fastq-dump -I --split-files ~/Downloads/lungCan_data_seurat/${sample}/${sample}.sra
# done



