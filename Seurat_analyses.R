library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./data") # pbmc data folder

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


pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# Normalizing layer: counts
# Performing log-normalization
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000) # with the default parameters

#'vst' method = relationship of log(variance) and log(mean) using loess
# Then standardizes the feature values using the observed mean and expected variance (given by the fitted line)
# Feature variance is calculated on the standardized values after clipping to a maximum
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

all.genes <- rownames(pbmc) # genes are rows
pbmc <- ScaleData(pbmc, features = all.genes)

# to remove unwanted sources of variation regress out heterogeneity
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc)) # runpca = pca for seurat object

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

# Runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique 
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
saveRDS(pbmc, file = "./data/output/pbmc_tutorial.rds")

#### Finding differentially expressed features (cluster biomarkers) ####

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)
# columns are p_val   avg_log2FC   pct.1   pct.2    p_val_adj

# finding all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)

# group pbmc.markers by cluster and keep rows that are avg_log2FC > 1
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# # A tibble: 7,024 × 7
# # Groups:   cluster [9]
# p_val  avg_log2FC  pct.1  pct.2  p_val_adj  cluster gene     


cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# plot raw counts as violin plots
VlnPlot(pbmc, features = c("NKG7", "PF4"), layer = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

# plot heatmap of top 10 genes
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
ggsave(filename = "./data/output/images/pbmc_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

