library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/") # pbmc data folder

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc", min.cells = 3, min.features = 200)
# pbmc

# The [[ operator can add columns to object metadata. This is a great place for QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# An object of class Seurat 
# 13714 features across 2638 samples within 1 assay 
# Active assay: RNA (13714 features, 0 variable features)
# 1 layer present: counts


pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# Normalizing layer: counts
# Performing log-normalization
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000) # with the default parameters

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Finding variable features for layer counts
# Calculating gene variances
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Calculating feature variances of standardized and clipped values
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
# to remove unwanted sources of variation
# we can regress out heterogeneity
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca") + NoLegend()

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

# Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
# To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
# This message will be shown once per session
# 10:12:00 UMAP embedding parameters a = 0.9922 b = 1.112
# 10:12:00 Read 2638 rows and found 10 numeric columns
# 10:12:00 Using Annoy for neighbor search, n_neighbors = 30
# 10:12:00 Building Annoy index with metric = cosine, n_trees = 50
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      10:12:00 Writing NN index file to temp file /var/folders/qw/jxd88n696kzf42ry24gryqwm0000gp/T//RtmpmKmrdB/filec03b9094c48
#    10:12:00 Searching Annoy index using 1 thread, search_k = 3000
#    10:12:01 Annoy recall = 100%
#    10:12:01 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
#    10:12:01 Initializing from normalized Laplacian + noise (using RSpectra)
#    10:12:01 Commencing optimization for 500 epochs, with 106310 positive edges
#    Using method 'umap'
#    0%   10   20   30   40   50   60   70   80   90   100%
#      [----|----|----|----|----|----|----|----|----|----|
#         **************************************************|
#         10:12:04 Optimization finished

DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "~/Downloads/filtered_gene_bc_matrices/output/pbmc_tutorial.rds")

#### Finding differentially expressed features (cluster biomarkers) ####

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)
# p_val avg_log2FC pct.1 pct.2    p_val_adj
# LTB  1.709675e-83   1.330256 0.982 0.647 2.344648e-79
# IL32 5.076510e-83   1.242930 0.947 0.471 6.961926e-79
# LDHB 2.467055e-68   1.044820 0.967 0.615 3.383320e-64
# CD3D 1.817480e-66   1.058609 0.920 0.438 2.492492e-62
# IL7R 8.698894e-61   1.389909 0.744 0.333 1.192966e-56

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# p_val avg_log2FC pct.1 pct.2     p_val_adj
# FCGR3A        5.972471e-204   6.795991 0.975 0.041 8.190647e-200
# IFITM3        5.671364e-195   6.201036 0.975 0.048 7.777708e-191
# CFD           2.389538e-193   6.081028 0.937 0.038 3.277012e-189
# CD68          1.800066e-189   5.472200 0.925 0.036 2.468611e-185
# RP11-290F20.3 6.852416e-189   6.390800 0.843 0.015 9.397404e-185

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# # A tibble: 7,024 × 7
# # Groups:   cluster [9]
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene     
# <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>    
#   1 5.32e-114       1.19 0.912 0.591 7.29e-110 0       LDHB     
# 2 1.31e- 83       2.35 0.439 0.11  1.79e- 79 0       CCR7     
# 3 2.61e- 78       1.06 0.85  0.403 3.58e- 74 0       CD3D     
# 4 5.89e- 55       1.03 0.731 0.398 8.07e- 51 0       CD3E     
# 5 3.91e- 50       2.11 0.338 0.104 5.36e- 46 0       LEF1     
# 6 2.53e- 47       1.23 0.624 0.36  3.47e- 43 0       NOSIP    
# 7 5.11e- 46       2.04 0.335 0.109 7.01e- 42 0       PRKCQ-AS1
# 8 5.49e- 43       1.51 0.438 0.186 7.52e- 39 0       PIK3IP1  
# 9 9.17e- 41       2.73 0.199 0.04  1.26e- 36 0       FHIT     
# 10 1.26e- 33       1.32 0.39  0.177 1.72e- 29 0       TCF7     
# # ℹ 7,014 more rows
# # ℹ Use `print(n = ...)` to see more rows

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

## Warning message:
##   The `slot` argument of `VlnPlot()` is deprecated as of Seurat 5.0.0.
## ℹ Please use the `layer` argument instead.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
# so use this instead:
VlnPlot(pbmc, features = c("NKG7", "PF4"), layer = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "~/Downloads/filtered_gene_bc_matrices/output/images/pbmc_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)





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


