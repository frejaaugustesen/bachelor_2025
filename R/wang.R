# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group
# QC of islet 34, NonDiabetic, male

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))
source(here::here("R/cell_types.R"))

set.seed(100)

wang_raw <- qread(here::here("data/seurat_objects/wang_sander_all_merged.qs"))
wang <- qread("/work/bachelor_2025/data/seurat_objects/wang_seurat_QC.qs")


# save QC  ----------------------------------------------------------------
# function to save QC as i go
qsave(wang, file = here::here("data/seurat_objects/wang_seurat_QC.qs"))
wang <- qread(here::here("data/seurat_objects/wang_seurat_QC.qs"))


wang <- qread("/work/bachelor_2025/data/seurat_objects/wang_sander_all_merged.qs")


# Integration with harmony ------------------------------------------------


# join count layers
wang[["RNA"]] <-SeuratObject::JoinLayers(wang[["RNA"]])

# split count layers by donor
wang[["RNA"]] <- base::split(wang[["RNA"]], f = wang$orig.ident)

# run standard anlaysis workflow
wang <- NormalizeData(wang)

wang <- FindVariableFeatures(wang, 
                                  selection.method = "vst", nfeatures = 2000)

# scale data
all.genes <- rownames(wang)
wang <- ScaleData(wang,
                       features = all.genes)


wang <- RunPCA(wang, features = VariableFeatures(object = wang))
 
ElbowPlot(wang) # 1:15

wang <- FindNeighbors(wang, dims = 1:15, reduction = "pca")
wang <- FindClusters(wang, resolution = 0.6)

wang <- RunUMAP(wang, dims = 1:15)

DimPlot(wang, reduction = "umap", label = TRUE)

# integration
wang <- IntegrateLayers(object = wang, 
                             method = HarmonyIntegration, 
                             orig.reduction = "pca", 
                             new.reduction = "HarmonyIntegration",
                             verbose = FALSE)

# re-join layers after integration
wang[["RNA"]] <- JoinLayers(wang[["RNA"]])

wang <- FindNeighbors(wang, reduction = "HarmonyIntegration", dims = 1:15)
wang <- FindClusters(wang, resolution = 0.6)

wang <- RunUMAP(wang, dims = 1:15, reduction = "HarmonyIntegration")


# New annotation ----------------------------------------------------------

## Clustering ----

DimPlot(wang, reduction = "umap", label = TRUE, group.by = "seurat_clusters")

## Dotplot ----
DotPlot(wang, 
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

DotPlot(wang, 
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
wang@meta.data <- wang@meta.data %>% 
  dplyr::mutate(wang_anno = dplyr::case_when(
    seurat_clusters %in% c(0, 3, 6, 8, 13) ~ "beta",
    seurat_clusters %in% c(1, 2, 4, 9) ~ "alpha",
    seurat_clusters %in% c(7) ~ "delta",
    seurat_clusters %in% c(10) ~ "gamma",
    seurat_clusters %in% c(5) ~"acinar",
    seurat_clusters %in% c(15) ~ "endothelial",
    seurat_clusters %in% c(12) ~ "stellate",
    seurat_clusters %in% c(14) ~ "immune",
    seurat_clusters %in% c(11) ~ "ductal",
  ))

DimPlot(wang, reduction = "umap", 
        label = TRUE, group.by = "wang_anno") + NoLegend()

# cell type distribution --------------------------------------------
table(wang@meta.data$wang_anno)
dark_palette <- colorRampPalette(brewer.pal(8, "Dark2"))(9)


# per sample
ggplot(wang@meta.data, aes(x = orig.ident, fill = wang_anno)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = dark_palette) +
  labs(
    title = "Cell type distribution per sample",
    x = "Sample",
    y = "Percentage"
  ) +
  theme(axis.text.x = element_text(angle = 45))



# agreement ---------------------------------------------------------------






