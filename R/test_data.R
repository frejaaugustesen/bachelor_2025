# test_data

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))

set.seed(100)

library(R.utils)

# load data ---------------------------------------------------------------

# unzip file
gunzip(here::here("test_data/pbmc3k_filtered_gene_bc_matrices.tar.gz"))

untar(here::here("test_data/pbmc3k_filtered_gene_bc_matrices.tar"))

data <- Read10X(here::here("filtered_gene_bc_matrices/hg19/"))

data_seurat <- CreateSeuratObject(counts = data, project = "pbmc3k", min.cells = 3, min.features = 200)

# project = det som objektet kommer til at stå som f.eks. donor id (kaldes orig.ident i meta data)


# data structure  ---------------------------------------------------------

str(data_seurat)

# QC and selecting cells --------------------------------------------------
## Violin plots, nFeatures, nCount_RNA og percent.mt ----

# tilføjer procent something - forstå lige hvad der bliver skrevet
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")

# Kigger på det procent der er blevet tilføjet
# samme måde at tilgå dataen - eventuelt brug str() til at se hvor den nye column ligger
data_seurat$percent.mt
data_seurat@meta.data$percent.mt

# plotter QC metrics
# Visualize QC metrics as a violin plot
VlnPlot(data_seurat, features = c("nFeature_RNA", 
                                  "nCount_RNA", "percent.mt"), ncol = 3)
# nFeature_RNA 
## hvor mange gener der udtrykkes i cellerne
## Den lange hale vil nok ud fra teori blive kategoriseret som "dubletter"
## Den korte hale/bunden vil være empty droplets, kun med ambient RNA

# nCount_RNA
## hvor mange unikke UMI der er??

# percent.mt - hvad er det?
## The percentage of reads that map to the mitochondrial genome
## We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, 
## which calculates the percentage of counts originating from a set of features
## We use the set of all genes starting with MT- as a set of mitochondrial genes

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

## feature-feature relations ----
plot1 <- FeatureScatter(data_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# hvad er det jeg kigger på?
head(data_seurat[["nFeature_RNA"]])
# måske mængden af gener som er fra RNA - altså ikke fra mt RNA

## filtering ----
data_seurat <- subset(data_seurat, 
               subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# filtrer celler, så celler som har færre end 200 og over 2500 gener fjernes,
# samt at celler med mere end 5% mitokondrielt RNA fjernes


# Normalization of data ---------------------------------------------------

data_seurat <- NormalizeData(data_seurat, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
# “LogNormalize” that normalizes the feature expression measurements for each cell 
# by the total expression, multiplies this by a scale factor (10,000 by default), 
# and log-transforms the result.

## view normalized data
data_seurat[["RNA"]]$data
## finder data "manuelt"
str(data_seurat)
data_seurat@assays$RNA@layers$data

# ovenstående i NormalizeData med method og scale.factor er default parametre
# havde givet det samme at skrive NormalizeData(data_seurat)

# Identification of highly variable features -----------------------------
data_seurat <- FindVariableFeatures(data_seurat,
                                    selection.method = "vst", 
                                    nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data_seurat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2 # Jeg kan ikke få lov til at plotte dem sammen
# virker fint hver for sig/individuelt
# er ikke helt sikker på at jeg forstår de her plots

# Scaling the data --------------------------------------------------------
## pre-processing prior to dimension reduction (PCA)
## Shifts the expression of each gene, so that the mean expression across cells is 0
## Scales the expression of each gene, so that the variance across cells is 1
## This step gives equal weight in downstream analyses, so that highly-expressed genes do 
## not dominate

all.genes <- rownames(data_seurat)
data_seurat <- ScaleData(data_seurat, features = all.genes)

head(data_seurat[["RNA"]]$scale.data)
# jeg forstår ikke hvad jeg kigger på?


# Perform linear dimensional reduction ------------------------------------

data_seurat <- RunPCA(data_seurat, features = VariableFeatures(object = data_seurat))
# Examine and visualize PCA results a few different ways
print(data_seurat[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(data_seurat, dims = 1:2, reduction = "pca")

DimPlot(data_seurat, reduction = "pca") + NoLegend()

DimHeatmap(data_seurat, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(data_seurat, dims = 1:15, cells = 500, balanced = TRUE) # maks 27 PC på en gang
str(data_seurat) 

data_seurat[["RunPCA.RNA"]]$npcs # der er 50 PC


# Determine the ‘dimensionality’ of the dataset ---------------------------

ElbowPlot(data_seurat)


# Cluster the cells -------------------------------------------------------


data_seurat <- FindNeighbors(data_seurat, dims = 1:10)
data_seurat <- FindClusters(data_seurat, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(data_seurat), 5)


# Run non-linear dimensional reduction (UMAP/tSNE) ------------------------

data_seurat <- RunUMAP(data_seurat, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(data_seurat, reduction = "umap") # mit vender omvendt i forhold til deres

# saving data
saveRDS(data_seurat, file = here::here("test_data/data_seurat.rds"))


# Finding differentially expressed features (cluster biomarkers) ----------

# find all markers of cluster 2
cluster2.markers <- FindMarkers(data_seurat, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(data_seurat, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
data_seurat.markers <- FindAllMarkers(data_seurat, only.pos = TRUE)
data_seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

saveRDS(data_seurat.markers, file = here::here("test_data/data_seurat.markers.rds"))
data_seurat.markers <- readRDS(here::here("test_data/data_seurat.markers.rds"))

cluster0.markers <- FindMarkers(data_seurat, ident.1 = 0, logfc.threshold = 0.25, 
                                test.use = "roc", only.pos = TRUE)

VlnPlot(data_seurat, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(data_seurat, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(data_seurat, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", 
                                      "FCGR3A", "LYZ", "PPBP", "CD8A"))

data_seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

DoHeatmap(data_seurat, features = top10$gene) + NoLegend()


# Assigning cell type identity to clusters --------------------------------

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(data_seurat)

data_seurat <- RenameIdents(data_seurat, new.cluster.ids)
DimPlot(data_seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(data_seurat, reduction = "umap", label = TRUE, label.size = 4.5) + 
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size = 10)))

ggsave(filename = here::here("test_data/final_plot.jpg"), height = 7, width = 12, plot = plot, quality = 50)

# histogramer

df <- as.data.frame(data_seurat@meta.data)
head(df)
ggplot(data = df, aes(x=percent.mt)) +
  geom_histogram(bins = 100)

