# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group

# first try

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))

set.seed(100)

# installing packages ----------------------------------------------------

## BiocManager 
install.packages("BiocManager")
library(BiocManager) # added to package_load

## scDblFinder - skal bruge hjælp
BiocManager::install("scDblFinder")
library(scDblFinder) # virker ikke??
BiocManager::install("plger/scDblFinder")

## Seurat
install.packages("Seurat")
library(Seurat) # added to package_load

## Jointly - skal bruge hjælp
install.packages("remotes")
remotes::install_github("madsen-lab/rJOINTLY")
library(rJOINTLY) # virker ikke??

## qs
install.packages("qs")
library(qs) # added to package_load

## XGBoost - understående er copy paste fra geeksforgeeks.org
# (XGBoost) Installing Packages 
install.packages("data.table") 
install.packages("dplyr") 
install.packages("ggplot2") 
install.packages("caret") 
install.packages("xgboost") 
install.packages("e1071") 
install.packages("cowplot") 

# (XGBoost) Loading packages 
library(data.table) # for reading and manipulation of data 
library(dplyr)     # for data manipulation and joining 
library(ggplot2) # for plotting 
library(caret)     # for modeling 
library(xgboost) # for building XGBoost model 
library(e1071)     # for skewness 
library(cowplot) # for combining multiple plots 

# loading sander seurat object -----------------------------------------------
sander <- qread("data/seurat_objects/wang_sander_merged.qs")

str(sander)
head(sander@assays$RNA@layers$counts)
sander[["RNA"]]$counts[1:10, 1:10] %>%
  as.matrix()

library(Seurat)
sander <- NormalizeData(sander)

head(sander@assays$RNA@layers$data)
sander[["RNA"]]$data[1:20, 1:10] %>%
  as.matrix()

str(sander)

sander[["RNA"]]


# motakis sample islet28 --------------------------------------------------

# using motakis Islet28
# folder for Islet 28: /data_raw/motakis/Islet28/Solo.out/GeneFull/raw

# loading data
islet28.mtx <- ReadMtx("data_raw/motakis/Islet28/Solo.out/GeneFull/raw/matrix.mtx", 
                   "data_raw/motakis/Islet28/Solo.out/GeneFull/raw/barcodes.tsv",
                   "data_raw/motakis/Islet28/Solo.out/GeneFull/raw/features.tsv")

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

islet28 <- CreateSeuratObject(counts.subset)

saveRDS(islet28, file = here::here("data/seurat_objects/motakis/islet28.rds"))

# Islet28 QC (motakis) -----------------------------------------------------¨

islet28[["percent.mt"]] <- PercentageFeatureSet(islet28, pattern = "^MT-")

## histograms ----
# creating dataframe of islet28 ti make histograms
islet28.df <- as.data.frame(islet28@meta.data)
head(islet28.df)

### nCount_RNA
ggplot(data = islet28.df, aes(x=nCount_RNA)) +
  geom_histogram(bins = 100)

### nFeatures_RNA
ggplot(data = islet28.df, aes(x=nFeature_RNA)) +
  geom_histogram(bins = 100)

### percent.mt
ggplot(data = islet28.df, aes(x=percent.mt)) +
  geom_histogram(bins = 100)

## subsetting data ----
# remember to asses histograms when choosing thresholds
islet28 <- subset(islet28, subset = nFeature_RNA > 900 &
                    nFeature_RNA < 7000 & percent.mt < 15 &
                    nCount_RNA > 9000, nCount_RNA < 40000)

## doublets ----

## marker genes ----

## clustering of cells ----




