# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group
# all of motakis UMAPs

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))

set.seed(100)


# load --------------------------------------------------------------------


# save objects ------------------------------------------------------------


# Merged before integration -----------------------------------------------

## seurat clusters ----
DimPlot(islet_all, reduction = "umap", label = TRUE, group.by = "seurat_clusters") +
  NoLegend()

## individual annotation ----
DimPlot(islet_all, reduction = "umap", label = TRUE, repel = TRUE,
        label.size = 3,
        group.by = "manual_anno") +
  NoLegend()

DimPlot(before, reduction = "umap", group.by = "orig.ident")



# Merged after integration ------------------------------------------------

## islet and seurat clusters ----
DimPlot(islet_all, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))
DimPlot(islet_all, reduction = "umap", group.by = "orig.ident")

DimPlot(wang, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(islet_all, reduction = "umap", group.by = "orig.ident")

## seurat clusters ----
DimPlot(islet_all, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  NoLegend()

## individual annotation ----
DimPlot(islet_all, reduction = "umap", repel = TRUE,
        label.size = 3,
        label = TRUE, 
        group.by = "manual_anno") +
  NoLegend()

## by disease ----
DimPlot(islet_all, reduction = "umap", group.by = "disease", pt.size = 0.1)

DimPlot(islet_all, split.by = "disease", reduction = "umap", group.by = "manual_anno",
        label = TRUE) + NoLegend()
DimPlot(islet_all, split.by = "disease", reduction = "umap", group.by = "disease")
DimPlot(islet_new, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  NoLegend()

DotPlot(islet_all, 
        features = azi_markers, group.by = "seurat_clusters"
)+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3.5, angle = 90))

# Doublets ----------------------------------------------------------------

## plotted doublets ----
DimPlot(islet_all, reduction = "umap", group.by = "multiplet")


# New annotation ----------------------------------------------------------
# annotation after integration and comparing with motakis annotation

DimPlot(islet_all, reduction = "umap", group.by = "merged_anno", label = TRUE) +
  NoLegend()

## new annotation ----

## motakis annotation ----

## agreement ----

## all three


# beta --------------------------------------------------------------------

motakis_beta <- subset(islet_all, subset = merged_anno == "beta")

DimPlot(motakis_beta, reduction = "umap", group.by = "merged_anno", label = TRUE) +
  NoLegend()

motakis_beta <- NormalizeData(motakis_beta)
motakis_beta <- FindVariableFeatures(motakis_beta,
                                     selection.method = "vst", nfeatures = 5000)
var_genes_motakis <- VariableFeatures(motakis_beta)
motakis_beta <- ScaleData(motakis_beta, features = var_genes_motakis)


# join count layers
motakis_beta[["RNA"]] <-SeuratObject::JoinLayers(motakis_beta[["RNA"]])

# split count layers by donor
motakis_beta[["RNA"]] <- base::split(motakis_beta[["RNA"]], f = motakis_beta$orig.ident)

# run standard anlaysis workflow
motakis_beta <- NormalizeData(motakis_beta)

motakis_beta <- FindVariableFeatures(motakis_beta, 
                                     selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(motakis_beta)

motakis_beta <- ScaleData(motakis_beta,
                          features = all.genes)


motakis_beta <- RunPCA(motakis_beta, features = VariableFeatures(object = motakis_beta))

ElbowPlot(motakis_beta) # 1:15