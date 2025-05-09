# working with all 12 samples 

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))
source(here::here("R/cell_types.R"))

set.seed(100)

# save QC ---------------------------
# newest version
qsave(islet_all, file = "/work/bachelor_2025/data/seurat_objects/motakis/islet_all_QC_new.qs")


# subset 
qsave(islet_all, file = "/work/bachelor_2025/data/seurat_objects/motakis/islet_all_subset.qs")
islet_all <- qread("/work/bachelor_2025/data/seurat_objects/motakis/islet_all_subset.qs")

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

# patcwork Dimplots ----
## ND ----
DimPlot57+DimPlot45+DimPlot38+DimPlot34

## preT2D ----
DimPlot56+DimPlot40+DimPlot63+DimPlot127

## T2D ----
DimPlot48+DimPlot52+DimPlot44+DimPlot128

# DotPlots --------------------------------------------------------------------

## 57 ----
DotPlot(
  islet57,
  features = azi_markers_short, group.by = "manual_anno"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("islet57 dotplot with Azimuth marker genes") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5),
    strip.text = element_text(angle = 45)
  )

## 45 ----
DotPlot(
  islet45,
  features = azi_markers_short, group.by = "manual_anno"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("islet45 dotplot with Azimuth marker genes") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5),
    strip.text = element_text(angle = 45)
  )

## 38 ----
DotPlot(
  islet38,
  features = azi_markers_short, group.by = "manual_anno"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("islet38 dotplot with Azimuth marker genes") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5),
    strip.text = element_text(angle = 45)
  )

## 34 ----
DotPlot(
  islet34,
  features = azi_markers_short, group.by = "manual_anno"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("islet34 dotplot with Azimuth marker genes") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5),
    strip.text = element_text(angle = 45)
  )

## 56 ----
DotPlot(
  islet56,
  features = azi_markers_short, group.by = "manual_anno"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("islet56 dotplot with Azimuth marker genes") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5),
    strip.text = element_text(angle = 45)
  )

## 40 ----
DotPlot(
  islet40,
  features = azi_markers_short, group.by = "manual_anno"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("islet40 dotplot with Azimuth marker genes") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5),
    strip.text = element_text(angle = 45)
  )

## 63 ----
DotPlot(
  islet63,
  features = azi_markers_short, group.by = "manual_anno"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("islet63 dotplot with Azimuth marker genes") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5),
    strip.text = element_text(angle = 45)
  )

## 127 ----
DotPlot(
  islet127,
  features = azi_markers_short, group.by = "manual_anno"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("islet127 dotplot with Azimuth marker genes") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5),
    strip.text = element_text(angle = 45)
  )

## 48 ----
DotPlot(
  islet48,
  features = azi_markers_short, group.by = "manual_anno"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("islet48 dotplot with Azimuth marker genes") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5),
    strip.text = element_text(angle = 45)
  )

## 52 ----
DotPlot(
  islet52,
  features = azi_markers_short, group.by = "manual_anno"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("islet52 dotplot with Azimuth marker genes") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5),
    strip.text = element_text(angle = 45)
  )

## 44 ----
DotPlot(
  islet44,
  features = azi_markers_short, group.by = "manual_anno"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("islet44 dotplot with Azimuth marker genes") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5),
    strip.text = element_text(angle = 45)
  )

## 128 ----
DotPlot(
  islet128,
  features = azi_markers_short, group.by = "manual_anno"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("islet128 dotplot with Azimuth marker genes") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5),
    strip.text = element_text(angle = 45)
  )


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
                         islet48, islet52, islet44, islet128),
                   add.cell.ids = c("islet57", "islet45", "islet38", "islet34",
                                    "islet56", "islet40", "islet63", "islet127",
                                    "islet48", "islet52", "islet44", "islet128"))

qsave(islet_all, file = "/work/bachelor_2025/data/seurat_objects/motakis/islet_all.qs")
islet_all <- qread("/work/bachelor_2025/data/seurat_objects/motakis/islet_all.qs")
before <- qread("/work/bachelor_2025/data/seurat_objects/motakis/islet_all.qs")

# integration with harmony ------------------------------------------------

# split the RNA measurements into two layers one for control cells, one for stimulated cells

islet_all
islet_all[["RNA"]]
str(islet_all)
islet_all@assays$RNA

# join count layers
islet_all[["RNA"]] <-SeuratObject::JoinLayers(islet_all[["RNA"]])

# split count layers by donor
islet_all[["RNA"]] <- base::split(islet_all[["RNA"]], f = islet_all$orig.ident)

# run standard anlaysis workflow
islet_all <- NormalizeData(islet_all)

islet_all <- FindVariableFeatures(islet_all, 
                                  selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(islet_all)

islet_all <- ScaleData(islet_all,
                       features = all.genes)


islet_all <- RunPCA(islet_all, features = VariableFeatures(object = islet_all))

ElbowPlot(islet_all) # 1:15

# islet_all <- FindNeighbors(islet_all, dims = 1:15, reduction = "pca")
# islet_all <- FindClusters(islet_all)

islet_all <- RunUMAP(islet_all, dims = 1:15)

DimPlot(islet_all, reduction = "umap", label = TRUE, group.by = "seurat_clusters") +
  NoLegend()

DimPlot(islet_all, reduction = "umap", label = TRUE, repel = TRUE,
        label.size = 3,
        group.by = "manual_anno") +
  NoLegend()

before <- islet_all

# integration
islet_all <- IntegrateLayers(object = islet_all, 
                             method = HarmonyIntegration, 
                             orig.reduction = "pca", 
                             new.reduction = "HarmonyIntegration",
                        verbose = FALSE)

# re-join layers after integration
islet_all[["RNA"]] <- JoinLayers(islet_all[["RNA"]])

islet_all <- FindNeighbors(islet_all, reduction = "HarmonyIntegration", dims = 1:15)
islet_all <- FindClusters(islet_all)

islet_all <- RunUMAP(islet_all, dims = 1:15, reduction = "HarmonyIntegration")

# Visualization
DimPlot(islet_all, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))
DimPlot(before, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))


DimPlot(islet_all, reduction = "umap", group.by = "seurat_clusters")


DimPlot(islet_all, reduction = "umap", repel = TRUE,
        label.size = 3,
        label = TRUE, 
        group.by = "manual_anno")

# tilføjer sygdom til meta data
meta <- read.csv(here::here("data_raw/motakis/motakis_meta.csv"), sep = ";") %>% 
  rename(orig.ident = "donor") %>% 
  mutate(orig.ident = tolower(orig.ident))
View(meta)

# add donor meta data to seurat object
islet_all@meta.data <-  islet_all@meta.data %>% 
  tibble::rownames_to_column("barcode") %>% 
  dplyr::left_join(y = meta, by = "orig.ident") %>% 
  tibble::column_to_rownames("barcode")

DimPlot(islet_all, reduction = "umap", group.by = "disease", pt.size = 0.1)

Idents(islet_all) <- "disease"

DimPlot(islet_all, split.by = "disease", reduction = "umap")

head(islet_all@meta.data)
str(islet_all@meta.data)

# isabell fikser problem når man har tilføjet meta data to gange
islet_all@meta.data <-  islet_all@meta.data %>% dplyr::select(-ends_with(".x"), -ends_with(".y"))


# doublets ----------------------------------------------------------------

## removing motakis data to avoid confusion

islet_all <- islet_all %>%
  Seurat::FindNeighbors(reduction = "pca",
                        dims = 1:15) %>%
  Seurat::FindClusters(resolution = 20,
                       algorithm = 1)

# Polyhormone detection - marker genes ------------------------------------

# Define marker genes to use for polyhormone detection
markers <- c("INS", "SST", "PPY", "GCG")

# find doublet clusters
# get average expression per cluster
avg <- Seurat::AverageExpression(islet_all, assay = "RNA",
                                 group.by = "RNA_snn_res.20")$RNA

# Scale gene expression column-wise (genes in columns)
avg.scaled <- t(scale(t(avg)))

# get polyhormone clusters

# Get scaled expression of canonical marker genes, and keep only clusters (columns) which have an
# average scaled expression above 0.6 (the expression is 0.6 standard deviation to the right of the mean on a bell curve (normal distribution).
# if the sum of these expression values within a cluster (column) is above 1 it means they expression
# more than 1 conical marker to a high degree, and thus could be doublets.
db_cluster <-
  names(which(colSums(avg.scaled[rownames(avg.scaled) %in% markers, ] > 0.6) > 1))

db_cluster_clean <- gsub("^g", "", db_cluster)

# save expression of marker genes
avg_scaled_df <- avg.scaled[rownames(avg.scaled) %in% markers, ] %>%
  as.data.frame() %>%
  dplyr::select(all_of(db_cluster))

# Plot doublets -----------------------------------------------------------
islet_all@meta.data <- islet_all@meta.data %>% 
  dplyr::mutate(multiplet = dplyr::case_when(RNA_snn_res.20 %in% db_cluster_clean ~ "multiplet",
                                             !RNA_snn_res.20 %in% db_cluster_clean ~ "singlet"))

## Doublet umap ----
DimPlot(islet_all, reduction = "umap", group.by = "multiplet")

FeaturePlot(islet_all, features = markers)

VlnPlot(islet_test, features = markers, group.by = "RNA_snn_res.20", 
        idents = db_cluster, pt.size = 0)

# remove doublets ----
islet_all <- subset(islet_all, subset = multiplet == "singlet")


# dotplot ----

islet_all <- FindClusters(islet_all, resolution = 0.5)

DotPlot(islet_all, 
                    features = azi_markers, group.by = "seurat_clusters"
                    )+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3.5, angle = 90))

# (Annotation) --------------------------------------------------------------

# cluster 0: alpha
# cluster 1: alpha
# cluster 2: alpha
# cluster 3 alpha
# cluster 4: beta
# cluster 5: alpha
# cluster 6: activated stellate
# cluster 7: acinar/ductal
# cluster 8: endothelial
# cluster 9: immune/ quiescent stellate

##  (annotering) ----

islet_all@meta.data <- islet_all@meta.data %>% 
  dplyr::mutate(merged_anno = dplyr::case_when(
  seurat_clusters %in% c(7) ~ "beta",
  seurat_clusters %in% c(8, 3, 2, 0) ~ "alpha",
  seurat_clusters %in% c(4) ~ "delta",
  seurat_clusters %in% c(11) ~"gamma",
  seurat_clusters %in% c(15) ~ "cycling",
  seurat_clusters %in% c(6) ~"acinar",
  seurat_clusters %in% c(14, 10) ~"endothelial",
  seurat_clusters %in% c(13) ~ "quiescent_stellate",
  seurat_clusters %in% c(5) ~"activated_stellate",
  seurat_clusters %in% c(12) ~"immune",
  seurat_clusters %in% c(9, 1) ~ "ductal",
), merged_anno = factor(merged_anno, levels = c("beta", "alpha","alpha_cycling",
                                            "delta", "gamma", "acinar", "endothelial",
                                            "quiescent_stellate", "activated_stellate",
                                            "immune", "ductal", "schwann")))

DimPlot(islet_test2, group.by = "seurat_clusters", reduction = "umap") 

DotPlot(
  islet_test2,
  features = azi_markers_short, group.by = "motakis_anno"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("Motakis annotation with Azimuth marker genes") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5)
  )

# replacing NA with "no annotation"
islet_all@meta.data <- islet_all@meta.data %>%
  mutate(merged_anno = ifelse(is.na(merged_anno), 
                              "no annotation", merged_anno))

# Feature Plot ---------------------------------------------------

FeaturePlot(islet_all, 
            features = c("INS", "GCG", "SST", "PPY"), pt.size = 0.1)


# comparing with motakis ----
install.packages("anndata")
library(anndata)

# loader h5ad objekt ind
motakis <- read_h5ad(here::here("data/Cellxgene Dataset.h5ad"))

#converter til seurat
motakis_s <- CreateSeuratObject(counts = t(as.matrix(motakis$X)), meta.data = motakis$obs)

motakis_meta <- motakis$obs
head(motakis_meta)

motakis_meta_sub <- motakis_meta %>% 
  tibble::rownames_to_column("barcode") %>% 
  mutate(Islet = tolower(Islet),
         Annotated_Clusters = tolower(Annotated_Clusters), 
         barcode = gsub(".*_","", barcode), 
         barcode = paste0(Islet, "_", barcode)) %>% 
  select(barcode, 
         motakis_anno = Annotated_Clusters, 
         motakis_nCount_RNA = nCount_RNA,
         motakis_nFeature_RNA = nFeature_RNA,
         motakis_percent.mt = percent.mt)

islet_df <- islet_all@meta.data

# motakis har valgt at lave to rækker for en celle hvis de ikke har kunne beslutte celletypen 
# derfor matcher jeg med første match (multiple = "first")
islet_df_2 <- islet_df %>% 
  tibble::rownames_to_column("barcode") %>% 
  left_join(y = motakis_meta_sub, by = "barcode", multiple = "first") %>% 
  mutate(agreement = case_when(manual_anno == motakis_anno ~ "yes", 
                               manual_anno != motakis_anno ~ "no")) %>% 
  tibble::column_to_rownames("barcode")
  

# regner procenter af yes/no i samlet df (islet_df_2)

islet_df_2 %>% 
  group_by(agreement) %>% 
  tally() %>% 
  mutate(perc_agree = round((n/sum(n))*100, 2))

islet_df_2 %>% 
  group_by(agreement) %>% 
  tally() %>% 
  mutate(perc_agree = round((n/sum(n))*100, 2),
         total = sum(n))

islet_all@meta.data <- islet_df_2

all.equal(colnames(islet_all), rownames(islet_df_2))

DimPlot(islet_all, reduction = "umap", label = TRUE,
        group.by = c("merged_anno", "motakis_anno", "agreement"),
        label.size = 2.5,
        repel = TRUE) & 
  NoLegend()


# 2 funktioner for det samme
table(islet_all@meta.data$agreement, islet_all@meta.data$manual_anno)

islet_all@meta.data %>% 
  group_by(agreement, manual_anno) %>% 
  tally()

#bruges ikke
islet_all@meta.data %>% 
  dplyr::filter(agreement == "no") %>% 
  pull("manual_anno") %>% 
  unique()
  

head(islet_all@meta.data)

# cell type distribution --------------------------------------------
table(islet_all@meta.data$merged_anno)

# per sample
ggplot(islet_all@meta.data, aes(x = orig.ident, fill = merged_anno)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Cell type distribution by sample",
    x = "Sample",
    y = "Percentage"
  ) +
  theme(axis.text.x = element_text(angle = 45))

# per disease
ggplot(islet_all@meta.data, aes(x = disease, fill = merged_anno)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Cell type distribution by disease state",
    x = "Disease state",
    y = "Percentage"
  )

# overall
ggplot(islet_all@meta.data, aes(x = "", fill = merged_anno)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Cell type distribution overall",
    x = "islet_all",
    y = "Percentage"
  )

head(islet_all@meta.data)

dark_palette <- colorRampPalette(brewer.pal(8, "Dark2"))(18)


# New annotation ----------------------------------------------------------

## Clustering ----
islet_all <- FindNeighbors(islet_all, dims = 1:15)
islet_all <- FindClusters(islet_all, resolution = 1)

islet_all <- RunUMAP(islet_all, dims = 1:15)

DimPlot(islet_all, reduction = "umap", label = TRUE, group.by = "seurat_clusters")

## Dotplot ----
DotPlot(islet_all, 
        features = azi_markers, group.by = "seurat_clusters"
)+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  ggtitle("Seurat clusters with Azimuth marker genes") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 3.5,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5),
    strip.text = element_text(angle = 45)
  )

## annotation ----
islet_all@meta.data <- islet_all@meta.data %>% 
  dplyr::mutate(merge_anno = dplyr::case_when(
  seurat_clusters %in% c(9) ~ "beta",
  seurat_clusters %in% c(0, 1, 2, 8, 10) ~ "alpha",
  seurat_clusters %in% c(7) ~ "delta",
  seurat_clusters %in% c(13) ~ "gamma",
  seurat_clusters %in% c(18) ~"cycling",
  seurat_clusters %in% c(4) ~"acinar",
  seurat_clusters %in% c(12) ~ "endothelial",
  seurat_clusters %in% c(16) ~ "quiescent_stellate",
  seurat_clusters %in% c(6, 17) ~ "activated_stellate",
  seurat_clusters %in% c(15) ~ "immune",
  seurat_clusters %in% c(3, 5, 11) ~ "ductal",
  seurat_clusters %in% c(14) ~ "schwann"
))

DimPlot(islet_all, reduction = "umap", 
        label = TRUE, group.by = "merge_anno") + NoLegend()

# nyt dimplot
DotPlot(islet_all, 
        features = azi_markers, group.by = "merge_anno"
)+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  ggtitle("Annotation with Azimuth marker genes") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 3.5,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5),
    strip.text = element_text(angle = 45)
  )

## agreement ----

# comparing with motakis 
reticulate::install_miniconda()

anndata::install_anndata()

library(anndata)

# loader h5ad objekt ind
motakis <- read_h5ad(here::here("data/Cellxgene Dataset.h5ad"))

#converter til seurat

motakis_meta <- motakis$obs
head(motakis_meta)

motakis_meta_sub <- motakis_meta %>% 
  tibble::rownames_to_column("barcode") %>% 
  mutate(Islet = tolower(Islet),
         Annotated_Clusters = tolower(Annotated_Clusters), 
         barcode = gsub(".*_","", barcode), 
         barcode = paste0(Islet, "_", barcode)) %>% 
  select(barcode, 
         motakis_anno = Annotated_Clusters, 
         motakis_nCount_RNA = nCount_RNA,
         motakis_nFeature_RNA = nFeature_RNA,
         motakis_percent.mt = percent.mt)

qsave(motakis_meta_sub, file = here::here("data/motakis_meta_sub.qs"))

# removing old motakis data
islet_df <- islet_all@meta.data %>% 
  select(-starts_with("motakis"))




# motakis har valgt at lave to rækker for en celle hvis de ikke har kunne beslutte celletypen 
# derfor matcher jeg med første match (multiple = "first")
islet_df <- islet_df %>% 
  tibble::rownames_to_column("barcode") %>% 
  left_join(y = motakis_meta_sub, by = "barcode", multiple = "first") %>% 
  mutate(agreement = case_when(merge_anno == motakis_anno ~ "yes", 
                               merge_anno != motakis_anno ~ "no")) %>% 
  tibble::column_to_rownames("barcode")

islet_df <- islet_df %>% 
  mutate(agreement = case_when(merge_anno == motakis_anno ~ "yes", 
                      merge_anno != motakis_anno ~ "no"))

# regner procenter af yes/no i samlet df (islet_df_2)

islet_df %>% 
  group_by(agreement) %>% 
  tally() %>% 
  mutate(perc_agree = round((n/sum(n))*100, 2))

agree_df <- islet_df %>% 
  group_by(merge_anno, agreement) %>% 
  tally() %>% 
  mutate(perc_agree = round((n/sum(n))*100, 2),
         total = sum(n))

islet_all@meta.data <- islet_df

all.equal(colnames(islet_all), rownames(islet_df_2))

DimPlot(islet_all, reduction = "umap", label = TRUE,
        group.by = c("merge_anno", "motakis_anno", "agreement"),
        label.size = 2.5,
        repel = TRUE) & 
  NoLegend()
 
 
# 2 funktioner for det samme
table(islet_all@meta.data$agreement, islet_all@meta.data$manual_anno)

# cell type distribution --------------------------------------------
table(islet_all@meta.data$manual_anno)

# per sample
ggplot(islet_all@meta.data, aes(x = orig.ident, fill = agreement)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Agreement by samples",
    x = "Sample",
    y = "Percentage"
  ) +
  theme(axis.text.x = element_text(angle = 45))


# per annotation
ggplot(islet_all@meta.data, aes(x = merge_anno, fill = agreement)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Agreement by annotation",
    x = "Disease state",
    y = "Percentage"
  ) +
  theme(axis.text.x = element_text(angle = 45))

# per motakis annotation
ggplot(islet_all@meta.data, aes(x = motakis_anno, fill = agreement)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Agreement by motakis annotation",
    x = "Disease state",
    y = "Percentage"
  ) +
  theme(axis.text.x = element_text(angle = 45))
