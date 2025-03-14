# working with all 12 samples 

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))

# loading in samples with QC ----
## ND ----
islet57 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet57_seurat_QC.rds"))
islet45 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet45_seurat_QC.rds"))
islet38 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet38_seurat_QC.rds"))
islet34 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet34_seurat_QC.rds"))

## preT2D ----
islet56 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet56_seurat_QC.rds"))
islet40 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet40_seurat_QC.rds"))
islet63 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet63_seurat_QC.rds"))
islet127 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet127_seurat_QC.rds"))

## T2D ----
islet48 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet48_seurat_QC.rds"))
islet52 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet52_seurat_QC.rds"))
islet44 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet44_seurat_QC.rds"))
islet128 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet128_seurat_QC.rds"))


# saving UMAPS with annotations ----

## ND ----
DimPlot57 <- DimPlot(islet57, reduction = "umap", label = TRUE, 
                     label.size = 3.5, pt.size = 0.5, repel = TRUE) + NoLegend()+
  ggtitle("islet57")

DimPlot45 <- DimPlot(islet45, reduction = "umap", label = TRUE, 
                     label.size = 3, pt.size = 0.5, repel = TRUE) + NoLegend()+
  ggtitle("islet45")

DimPlot38 <- DimPlot(islet38, reduction = "umap", label = TRUE, 
                     label.size = 3, pt.size = 0.5, repel = TRUE) + NoLegend()+
  ggtitle("islet38")

DimPlot34 <- DimPlot(islet34, reduction = "umap", label = TRUE, 
                     label.size = 3, pt.size = 0.5, repel = TRUE) + NoLegend() +
  ggtitle("islet34")

## preT2D ----

DimPlot56 <- DimPlot(islet56, reduction = "umap", label = TRUE, 
                     label.size = 3, pt.size = 0.5, repel = TRUE) + NoLegend()+
  ggtitle("islet56")

DimPlot40 <- DimPlot(islet40, reduction = "umap", label = TRUE, 
                     label.size = 3, pt.size = 0.5, repel = TRUE) + NoLegend()+
  ggtitle("islet40")

DimPlot63 <- DimPlot(islet63, reduction = "umap", label = TRUE, 
                     label.size = 3, pt.size = 0.5, repel = TRUE) + NoLegend()+
  ggtitle("islet63")

DimPlot127 <- DimPlot(islet127, reduction = "umap", label = TRUE, 
                     label.size = 3, pt.size = 0.5, repel = TRUE) + NoLegend()+
  ggtitle("islet127")

## T2D ----

DimPlot48 <- DimPlot(islet48, reduction = "umap", label = TRUE, 
                     label.size = 3, pt.size = 0.5, repel = TRUE) + NoLegend()+
  ggtitle("islet48")

DimPlot52 <- DimPlot(islet52, reduction = "umap", label = TRUE, 
                     label.size = 3, pt.size = 0.5, repel = TRUE) + NoLegend()+
  ggtitle("islet52")

DimPlot44 <- DimPlot(islet44, reduction = "umap", label = TRUE, 
                     label.size = 3, pt.size = 0.5, repel = TRUE) + NoLegend()+
  ggtitle("islet44")

DimPlot128 <- DimPlot(islet128, reduction = "umap", label = TRUE, 
                      label.size = 3, pt.size = 0.5, repel = TRUE) + NoLegend()+
  ggtitle("islet128")

# patcwork plots ----
## ND ----
DimPlot57+DimPlot45+DimPlot38+DimPlot34

## preT2D ----
DimPlot56+DimPlot40+DimPlot63+DimPlot127

## T2D ----
DimPlot48+DimPlot52+DimPlot44+DimPlot128


# Merged seurat object ----------------------------------------------------

islet_ND <- merge(islet57, 
                  y= c(islet45, islet38, islet34),
                  add.cell.ids = c("islet57", "islet45", "islet38", "islet34"))

islet_pre <- merge(islet56, 
                   y= c(islet40, islet63, islet127),
                   add.cell.ids = c("islet56", "islet40", "islet63", "islet127"))

islet_T2D <- merge(islet48, 
                   y= c(islet52, islet44, islet128),
                   add.cell.ids = c("islet48", "islet52", "islet44", "islet128"))

islet_all <- merge(islet57,
                   y= c(islet45, islet38, islet34, 
                         islet56, islet40, islet63, islet127,
                         islet57, islet52, islet44, islet128),
                   add.cell.ids = c("islet57", "islet45", "islet38", "islet34",
                                    "islet56", "islet40", "islet63", "islet127",
                                    "islet48", "islet52", "islet44", "islet128"))
saveRDS(islet_all, file = "/work/bachelor_2025/data/seurat_objects/motakis/islet_all.rds")


# integration with harmony ------------------------------------------------

# split the RNA measurements into two layers one for control cells, one for stimulated cells

islet_all
islet_all[["RNA"]]
str(islet_all)
islet_all@assays$RNA

islet_all[["RNA"]] <- split(islet_all[["RNA"]], f = islet_all$orig.ident)

# run standard anlaysis workflow
islet_all <- NormalizeData(islet_all)
islet_all <- FindVariableFeatures(islet_all)
islet_all <- ScaleData(islet_all)
islet_all <- RunPCA(islet_all, features = VariableFeatures(object = islet_all))

ElbowPlot(islet_all) # 1:15

islet_all <- FindNeighbors(islet_all, dims = 1:15, reduction = "pca")
islet_all <- FindClusters(islet_all, resolution = 0.6)

islet_all <- RunUMAP(islet_all, dims = 1:15)

DimPlot(islet_all, reduction = "umap", label = TRUE)

DimPlot(islet_all, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))

# integration
islet_all <- IntegrateLayers(object = islet_all, 
                             method = HarmonyIntegration, 
                             orig.reduction = "pca", 
                             new.reduction = "HarmonyIntegration",
                        verbose = FALSE)

# re-join layers after integration
islet_all[["RNA"]] <- JoinLayers(islet_all[["RNA"]])

islet_all <- FindNeighbors(islet_all, reduction = "HarmonyIntegration", dims = 1:15)
islet_all <- FindClusters(islet_all, resolution = 0.6)

islet_all <- RunUMAP(islet_all, dims = 1:15, reduction = "HarmonyIntegration")

# Visualization
DimPlot(islet_all, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))

DimPlot(islet_all, reduction = "umap", label = TRUE)

