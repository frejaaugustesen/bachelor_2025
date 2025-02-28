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

# find empty droplets
remotes::install_github("madsen-lab/valiDrops")
library(valiDrops)

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
# creating dataframe of islet28 ti make histograms
islet28.df <- as.data.frame(islet28@meta.data)
head(islet28.df)

### nCount_RNA
ggplot(data = islet28.df, aes(x=nCount_RNA)) +
  geom_histogram(bins = 200) + 
  geom_vline(xintercept = 1500)


### nFeatures_RNA
ggplot(data = islet28.df, aes(x=nFeature_RNA) +
         geom_histogram(bins = 100)))
  

### percent.mt
ggplot(data = islet28.df, aes(x=percent.mt)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 15)


## subsetting data ----
# remember to asses histograms when choosing thresholds
islet28 <- subset(islet28, subset = nFeature_RNA >= 1000 & 
                    percent.mt < 15 &
                    nCount_RNA >= 1500 & nCount_RNA <= 40000)

## kompleksitet af celler (nFeature/nCount) ----
islet28@meta.data$log10complexity <- log10(islet28@meta.data$nFeature_RNA)/
  log10(islet28@meta.data$nCount_RNA)

hist(islet28@meta.data$log10complexity, )

ggplot(data = islet28@meta.data, aes(x=log10complexity)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 0.80, col = "red")

## doublets (udskudt) ----

## marker genes ----

## clustering of cells ----


#  pancreas azimuth -------------------------------------------------------

panc_data <-readRDS("test_data/enge.rds")
str(panc_data)

panc_df <- as.data.frame(panc_data@meta.data)

ggplot(data = panc_df, aes(x=nCount_RNA)) +
  geom_histogram(bins = 100)

ggplot(data = panc_df, aes(x=nFeature_RNA)) +
  geom_histogram(bins = 100)

head(panc_data@assays$RNA$data, 10)

# using guide

panc_data <- FindVariableFeatures(panc_data, 
                                  selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(panc_data), 10)

plot1 <- VariableFeaturePlot(panc_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# scaling data
all.genes <- rownames(panc_data)
panc_data <- ScaleData(panc_data, features = all.genes)

# pca
panc_data <- RunPCA(panc_data, features = VariableFeatures(object = panc_data))

ElbowPlot(panc_data)

# clustering
panc_data <- FindNeighbors(panc_data, dims = 1:13)
panc_data <- FindClusters(panc_data, resolution = 0.5)

panc_data <- RunUMAP(panc_data, dims = 1:13)

DimPlot(panc_data, reduction = "umap")

