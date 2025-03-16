# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group
# QC of islet 34, NonDiabetic, male

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))
source(here::here("R/cell_types.R"))

set.seed(100)

islet34_raw <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet34_seurat.rds"))
islet34 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet34_seurat_QC.rds"))


# save QC  ----------------------------------------------------------------
# function to save QC as i go
saveRDS(islet34, file = here::here("data/seurat_objects/motakis/selected_samples/islet34_seurat_QC.rds"))




# load data ---------------------------------------------------------------

# loader data ind som matrix (rækker = gener, kolonner = celler)
mtx_gene <- Seurat::ReadMtx(mtx = "data_raw/motakis/Islet34/Solo.out/Gene/raw/matrix.mtx",
                            cells = "data_raw/motakis/Islet34/Solo.out/Gene/raw/barcodes.tsv",
                            features = "data_raw/motakis/Islet34/Solo.out/Gene/raw/features.tsv")

mtx_genefull <- Seurat::ReadMtx(mtx = "data_raw/motakis/Islet34/Solo.out/GeneFull/raw/matrix.mtx",
                                cells = "data_raw/motakis/Islet34/Solo.out/GeneFull/raw/barcodes.tsv",
                                features = "data_raw/motakis/Islet34/Solo.out/GeneFull/raw/features.tsv")

# threshold har en df ("ranks") med 3 kolonner, hvor den første er barcode, næste er count 
# (faldende) og sidste er rank (stigende), samt en værdi der hedder "lower.threshold" (647)
threshold <- valiDrops::rank_barcodes(mtx_gene) 

# rank.pass er en liste over de barcodes som er over/= "lower.threshold" værdien
# man beholder rownames fra df "ranks", hvis værdien i "counts"-kolonne er >= lower threshold
rank.pass <- BiocGenerics::rownames(threshold$ranks
                                    [ threshold$ranks$counts >= threshold$lower.threshold, ])

# tager de kolonne-navne (celler) fra mtx gene som også er tilstede i rank.pass.
# subsettet matrix med kun de celler som er over threshold
mtx_gene_sub <- mtx_gene[, colnames(mtx_gene) %in% rank.pass]

# subsetter matrix med de navne som også er til stede i mtx_gene_sub
# ville det her give det samme at subsette med rank.pass?
# nu indeholder de to sub matrix ihvertfald kun de samme celler
mtx_genefull_sub <- mtx_genefull[, colnames(mtx_genefull) %in% colnames(mtx_gene_sub)]

# calcuclate exon-intron / exon ratio
# tester at alle celler er de samme i begge sub matrix - havde givet false hvis de ikke var
all.equal(colnames(mtx_gene_sub), colnames(mtx_genefull_sub))

# colSums(mtx_gene_sub) = samlet antal gener udtrykt i hver celle
# samlede antal gener udtrykt i hver celle fra gene og genefull divideret og lavet til en df
# dette er exon / exon-intron ratioen
contrast <- colSums(mtx_gene_sub) / colSums(mtx_genefull_sub) %>% as.data.frame()

# kollonen med ratioen navngives
colnames(contrast) <- "exon_exonintron_ratio"

# laver seurat objekt med sub matrix
islet34 <- Seurat::CreateSeuratObject(counts = mtx_gene_sub,
                                      assay = "RNA",
                                      project = "islet34")

# tjekker at der er de samme celler i contrast dataframe og i seurat objektet
all.equal(rownames(contrast), rownames(islet34@meta.data))

# tilføjer ratioen til meta data
islet34@meta.data$exon_exonintron <- contrast$exon_exonintron

# gemmer seuratobjekt
saveRDS(islet34, file = here::here("data/seurat_objects/motakis/selected_samples/islet34_seurat.rds"))


# QC ----------------------------------------------------------------------

# adding percentage of mitochondrial RNA ----
islet34[["percent.mt"]] <- PercentageFeatureSet(islet34, pattern = "^MT-")

# kompleksitet af celler (nFeature/nCount) ----
islet34@meta.data$log10complexity <- log10(islet34@meta.data$nFeature_RNA)/
  log10(islet34@meta.data$nCount_RNA)

## histograms ----

### nCount_RNA
ggplot(data = islet34@meta.data, aes(x=nCount_RNA)) +
  geom_histogram(bins = 200) + 
  geom_vline(xintercept = 7000) +
  geom_vline(xintercept = 40000)


### nFeatures_RNA
ggplot(data = islet34@meta.data, aes(x=nFeature_RNA)) +
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = 1900)


### percent.mt
ggplot(data = islet34@meta.data, aes(x=percent.mt)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 15)

### exon intron ratio
ggplot(data = islet34@meta.data, aes(x=exon_exonintron)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 1)

### complexity
ggplot(data = islet34@meta.data, aes(x=log10complexity)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 0.80)

## subsetting data ----
# remember to asses histograms when choosing thresholds
islet34 <- subset(islet34, subset = nFeature_RNA >= 1900 & 
                    percent.mt < 15 &
                    nCount_RNA >= 7000 & nCount_RNA <= 40000 & 
                    exon_exonintron <= 1 & log10complexity >= 0.80)

# saving QC as seurat object - to get the previous histograms, load raw seurat
saveRDS(islet34, file = here::here("data/seurat_objects/motakis/selected_samples/islet34_seurat_QC.rds"))

## normalizing data ----
islet34 <- NormalizeData(islet34)

## variable features ----
islet34 <- FindVariableFeatures(islet34,
                                selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(islet34), 10)

# plot variable features with and without labels
varfeat_plot_34 <- VariableFeaturePlot(islet34)
varfeat_plot_34_2 <- LabelPoints(plot = varfeat_plot_34,
                                 points = top10, repel = TRUE)

## scaling data and PCA ----
all.genes <- rownames(islet34)
islet34 <- ScaleData(islet34, features = all.genes)

# PCA
islet34 <- RunPCA(islet34, features = VariableFeatures(object = islet34))

# assesing important PCs
ElbowPlot(islet34) # 1:15


# clustering  -------------------------------------------------------------

## marker genes ----
islet34 <- FindNeighbors(islet34, dims = 1:15)
islet34 <- FindClusters(islet34, resolution = 0.6)

## clustering of cells ----
islet34 <- RunUMAP(islet34, dims = 1:15)

DimPlot(islet34, reduction = "umap", label = TRUE)

# Dotplots ----------------------------------------------------------------


DotPlot1 <- DotPlot(islet34, 
                    features = list("beta"=beta, "alpha" = alpha, 
                                    "delta" = delta, "gamma" = gamma
                    ))+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3.5))

DotPlot2 <- DotPlot(islet34, 
                    features = list("acinar" = acinar, "ductal" = ductal, 
                                    "cycling" = cycling, "immune" = immune
                    ))+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3.5))

DotPlot3 <- DotPlot(islet34, 
                    features = list("activated_stellate" = activated_stellate,
                                    "endothelial" = endothelial, "epsilon" = epsilon,
                                    "quiescent_stellate" = quiescent_stellate
                    ))+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3.5))


# Annotation --------------------------------------------------------------

# cluster 0: alpha
# cluster 1: beta
# cluster 2: ductal
# cluster 3: ductal
# cluster 4: acinar
# cluster 5: delta
# cluster 6: activated stellate
# cluster 7: beta/alpha
# cluster 8: gamma
# cluster 9: immune
# cluster 10:cycling
# cluster 11: quiescent stellate/endothelial

## klade til annotering ----

islet34@meta.data <- islet34@meta.data %>% dplyr::mutate(manual_anno = dplyr::case_when(
  seurat_clusters %in% c(0) ~ "alpha",
  seurat_clusters %in% c(7) ~ "beta/alpha",
  seurat_clusters %in% c(6) ~ "activated_stellate",
  seurat_clusters %in% c(4) ~ "acinar",
  seurat_clusters %in% c(11) ~"quiescent_stellate/endothelial",
  seurat_clusters %in% c(9) ~"immune",
  seurat_clusters %in% c(10) ~ "cycling",
  seurat_clusters %in% c(2, 3) ~ "ductal",
  seurat_clusters %in% c(8) ~ "gamma",
  seurat_clusters %in% c(5) ~ "delta",
  seurat_clusters %in% c(1) ~ "beta",
))

DimPlot(islet34, group.by = "manual_anno", 
        reduction = "umap", label = TRUE, repel = TRUE) +
  NoLegend()

DimPlot(islet34, group.by = "seurat_clusters", 
        reduction = "umap", label = TRUE, repel = TRUE) +
  NoLegend()



# small dotplots ----------------------------------------------------------
alpha1 <- c("GCG", "TTR", "PEMT")
beta1 <- c("IAPP", "INS", "NPTX2")
ductal1 <- c("SPP1", "SERPINA1", "CFTR")
acinar1 <- c("REG1A", "PRSS1", "CTRB2")
delta1 <- c("SST", "RBP4", "SEC11C")
activated_stellate1 <- c("COL1A1", "COL1A2", "BGN")
gamma1 <- c("PPY", "ETV1", "MEIS2")
immune1 <- c("ACP5", "APOE", "C1QB")
cycling1 <- c("UBE2C", "BIRC5", "CDKN3")
quiescent_stellate1 <- c("ESAM", "IGFBP4", "CSRP2")
endothelial1 <- c("PLVAP", "RGCC", "PECAM1")

DotPlot(islet38, features = list("beta"=beta1, "alpha" = alpha1, "delta" = delta1,
                                 "acinar" = acinar1, "ductal" = ductal1,
                                 "activated_stellate" = activated_stellate1,
                                 "endothelial" = endothelial1, "immune" = immune1, 
                                 "gamma" = gamma1, "cycling" = cycling1,
                                 "quiescent_stellate" = quiescent_stellate1
))+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 4, 
                                   angle = 90, vjust = 0.5, hjust = 0.5))


# Heatmap -----------------------------------------------------------------

DoHeatmap(islet34, features = c(beta, alpha, delta, gamma, epsilon, cycling, ductal, 
                                endothelial, immune, quiescent_stellate,
                                schwann, activated_stellate, acinar), size = 2) +
  theme(text = element_text(size = 6))


# Variable Feature Plot ---------------------------------------------------

FeaturePlot(islet34, features = c("INS", "GCG", "SST", "PPY"), pt.size = 1)

