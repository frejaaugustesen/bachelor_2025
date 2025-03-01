# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group

# first try

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))

set.seed(100)

# islet28 load (motakis) --------------------------------------------------

# using motakis Islet28
# folder for Islet 28: /data_raw/motakis/Islet28/Solo.out/GeneFull/raw

# loading data
islet28.mtx <- ReadMtx("data_raw/motakis/Islet28/Solo.out/Gene/raw/matrix.mtx", 
                   "data_raw/motakis/Islet28/Solo.out/Gene/raw/barcodes.tsv",
                   "data_raw/motakis/Islet28/Solo.out/Gene/raw/features.tsv")

# creating plot and threshold
threshold <- valiDrops::rank_barcodes(islet28.mtx, type = "UMI") # tager virkelig lang tid!

# finding cells that pass threshold
rank.pass <- BiocGenerics::rownames(threshold$ranks[ threshold$ranks$counts >= threshold$lower.threshold, ])

# subsetting matrix
counts.subset <- islet28.mtx[, colnames(islet28.mtx) %in% rank.pass]

# creating seurat object

islet28 <- CreateSeuratObject(counts.subset, project = "islet28") #HUSK

saveRDS(islet28, file = here::here("data/seurat_objects/motakis/islet28.rds"))
islet28 <- readRDS(here::here("data/seurat_objects/motakis/islet28.rds"))

# Islet28 QC (motakis) -----------------------------------------------------

islet28[["percent.mt"]] <- PercentageFeatureSet(islet28, pattern = "^MT-")

## histograms ----

### nCount_RNA
ggplot(data = islet28@meta.data, aes(x=nCount_RNA)) +
  geom_histogram(bins = 200) + 
  geom_vline(xintercept = 1200)


### nFeatures_RNA
ggplot(data = islet28@meta.data, aes(x=nFeature_RNA)) +
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = 800)
  

### percent.mt
ggplot(data = islet28@meta.data, aes(x=percent.mt)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 15)


## subsetting data ----
# remember to asses histograms when choosing thresholds
islet28 <- subset(islet28, subset = nFeature_RNA >= 800 & 
                    percent.mt < 15 &
                    nCount_RNA >= 1200 & nCount_RNA <= 40000)

## kompleksitet af celler (nFeature/nCount) ----
islet28@meta.data$log10complexity <- log10(islet28@meta.data$nFeature_RNA)/
  log10(islet28@meta.data$nCount_RNA)

ggplot(data = islet28@meta.data, aes(x=log10complexity)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 0.80, col = "red")

## vionlinplots overview ----
VlnPlot(islet28, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## normalizing data ----
islet28 <- NormalizeData(islet28)

## variable features ----
islet28 <- FindVariableFeatures(islet28, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(islet28), 10)

# plot variable features with and without labels
varfeat_plot_28 <- VariableFeaturePlot(islet28)
varfeat_plot_28_2 <- LabelPoints(plot = varfeat_plot_28, points = top10, repel = TRUE)

## scaling data and PCA ----
all.genes <- rownames(islet28)
islet28 <- ScaleData(islet28, features = all.genes)

# PCA
islet28 <- RunPCA(islet28, features = VariableFeatures(object = islet28))

# assesing important PCs
ElbowPlot(islet28) # 1:16

## marker genes ----
islet28 <- FindNeighbors(islet28, dims = 1:16)
islet28 <- FindClusters(islet28, resolution = 0.6)

head(Idents(islet28), 5)

## clustering of cells ----
islet28 <- RunUMAP(islet28, dims = 1:16)

DimPlot(islet28, reduction = "umap")

# find markers
islet28.markers <- FindAllMarkers(islet28, only.pos = TRUE)

saveRDS(islet28.markers, file = here::here("data/rds_files/islet28/islet28.markers.rds"))
saveRDS(islet28, file = here::here("data/seurat_objects/motakis/islet28/islet28_QC.rds"))

islet28.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

VlnPlot(islet28, features = c("CYSTM1", "CRYBA2"))

# FeaturePlots

FeaturePlot(islet28, features = "RPL21")

FeaturePlot(islet28, features = beta)

# markers for each cluster (tror at cluster 1, 4 og 6 er beta)
# cluster nul anderledes 1 og 2
cluster0.markers <- FindMarkers(islet28, ident.1 = 0)
head(cluster0.markers, n = 5) # måske PPP1R1A (alfa)

# beta
beta_markers_in_clusters <- islet28.markers[islet28.markers$gene 
                                            %in% beta, ]
table(beta_markers_in_clusters$cluster) # 1, 11 og 12 (6 og 4)

# aplha
alpha_markers_in_clusters <- islet28.markers[islet28.markers$gene 
                                            %in% alpha, ]
table(alpha_markers_in_clusters$cluster) # 0 (3)

# delta
delta_markers_in_clusters <- islet28.markers[islet28.markers$gene 
                                             %in% delta, ]
table(delta_markers_in_clusters$cluster) # 5 og 12

# gamma
gamma_markers_in_clusters <- islet28.markers[islet28.markers$gene 
                                             %in% gamma, ]
table(gamma_markers_in_clusters$cluster) # 5


# random sample selection -------------------------------------------------

meta <- read.csv(here::here("data_raw/motakis/motakis_meta.csv"), sep = ";")
View(meta)

table(meta$disease)

samples <- meta %>%
  group_by(disease) %>% # grupperer data
  sample_n(3)  %>% # udvælger tilfældigt 3 samples
  ungroup() 

print(samples %>% select(donor, disease))
