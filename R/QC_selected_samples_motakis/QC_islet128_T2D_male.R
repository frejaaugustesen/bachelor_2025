# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group
# QC of islet 128, T2Diabetic, male

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))
source(here::here("R/cell_types.R"))

set.seed(100)

islet128_raw <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet128_seurat.rds"))
islet128 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet128_seurat_QC.rds"))


# save QC  ----------------------------------------------------------------
# function to save QC as i go
saveRDS(islet128, file = here::here("data/seurat_objects/motakis/selected_samples/islet128_seurat_QC.rds"))


# load data ---------------------------------------------------------------

# loader data ind som matrix (rækker = gener, kolonner = celler)
mtx_gene <- Seurat::ReadMtx(mtx = "data_raw/motakis/Islet128/Solo.out/Gene/raw/matrix.mtx",
                            cells = "data_raw/motakis/Islet128/Solo.out/Gene/raw/barcodes.tsv",
                            features = "data_raw/motakis/Islet128/Solo.out/Gene/raw/features.tsv")

mtx_genefull <- Seurat::ReadMtx(mtx = "data_raw/motakis/Islet128/Solo.out/GeneFull/raw/matrix.mtx",
                                cells = "data_raw/motakis/Islet128/Solo.out/GeneFull/raw/barcodes.tsv",
                                features = "data_raw/motakis/Islet128/Solo.out/GeneFull/raw/features.tsv")

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
islet128 <- Seurat::CreateSeuratObject(counts = mtx_gene_sub,
                                      assay = "RNA",
                                      project = "islet128")

# tjekker at der er de samme celler i contrast dataframe og i seurat objektet
all.equal(rownames(contrast), rownames(islet128@meta.data))

# tilføjer ratioen til meta data
islet128@meta.data$exon_exonintron <- contrast$exon_exonintron

# gemmer seuratobjekt
saveRDS(islet128, file = here::here("data/seurat_objects/motakis/selected_samples/islet128_seurat.rds"))


# QC ----------------------------------------------------------------------

# adding percentage of mitochondrial RNA ----
islet128[["percent.mt"]] <- PercentageFeatureSet(islet128, pattern = "^MT-")

# kompleksitet af celler (nFeature/nCount) ----
islet128@meta.data$log10complexity <- log10(islet128@meta.data$nFeature_RNA)/
  log10(islet128@meta.data$nCount_RNA)

## histograms ----

### nCount_RNA
ggplot(data = islet128@meta.data, aes(x=nCount_RNA)) +
  geom_histogram(bins = 200) + 
  geom_vline(xintercept = 4000) +
  geom_vline(xintercept = 50000)


### nFeatures_RNA
ggplot(data = islet128@meta.data, aes(x=nFeature_RNA)) +
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = 2500)


### percent.mt
ggplot(data = islet128@meta.data, aes(x=percent.mt)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 25)

### exon intron ratio
ggplot(data = islet128@meta.data, aes(x=exon_exonintron)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 1)

### complexity
ggplot(data = islet128@meta.data, aes(x=log10complexity)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 0.80)

## subsetting data ----
# remember to asses histograms when choosing thresholds
islet128 <- subset(islet128, subset = nFeature_RNA >= 2500 & 
                    percent.mt < 25 &
                    nCount_RNA >= 4000 & nCount_RNA <= 50000 & 
                    exon_exonintron <= 1 & log10complexity >= 0.80)

# saving QC as seurat object - to get the previous histograms, load raw seurat
saveRDS(islet128, file = here::here("data/seurat_objects/motakis/selected_samples/islet128_seurat_QC.rds"))

## normalizing data ----
islet128 <- NormalizeData(islet128)

## variable features ----
islet128 <- FindVariableFeatures(islet128,
                                selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(islet128), 10)

# plot variable features with and without labels
varfeat_plot_128 <- VariableFeaturePlot(islet128)
varfeat_plot_128_2 <- LabelPoints(plot = varfeat_plot_128,
                                 points = top10, repel = TRUE)

## scaling data and PCA ----
all.genes <- rownames(islet128)
islet128 <- ScaleData(islet128, features = all.genes)

# PCA
islet128 <- RunPCA(islet128, features = VariableFeatures(object = islet128))

# assesing important PCs
ElbowPlot(islet128) # 1:13

# clustering  -------------------------------------------------------------

## marker genes ----
islet128 <- FindNeighbors(islet128, dims = 1:13)
islet128 <- FindClusters(islet128, resolution = 0.6)

## clustering of cells ----
islet128 <- RunUMAP(islet128, dims = 1:13)

DimPlot(islet128, reduction = "umap", label = TRUE)

# Dotplots ----------------------------------------------------------------

DotPlot1 <- DotPlot(islet128, 
                    features = list("beta"=beta, "alpha" = alpha, 
                                    "delta" = delta, "gamma" = gamma
                    ))+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3.5))

DotPlot2 <- DotPlot(islet128, 
                    features = list("acinar" = acinar, "ductal" = ductal, 
                                    "cycling" = cycling, "immune" = immune
                    ))+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3.5))

DotPlot3 <- DotPlot(islet128, 
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
# cluster 1: ductal
# cluster 2: gamma
# cluster 3: activated_stellate
# cluster 4: alpha
# cluster 5: beta
# cluster 6: delta
# cluster 7: acinar
# cluster 8: endothelial
# cluster 9: quiescent stellate
# cluster 10: immune


## klade til annotering ----

islet128@meta.data <- islet128@meta.data %>% dplyr::mutate(manual_anno = dplyr::case_when(
  seurat_clusters %in% c(0, 4) ~ "alpha",
  seurat_clusters %in% c(5) ~ "beta",
  seurat_clusters %in% c(3) ~ "activated_stellate",
  seurat_clusters %in% c(7) ~ "acinar",
  seurat_clusters %in% c(1) ~ "ductal",
  seurat_clusters %in% c(6) ~ "delta",
  seurat_clusters %in% c(8) ~ "endothelial",
  seurat_clusters %in% c(10) ~ "immune",
  seurat_clusters %in% c(9) ~ "quiescent_stellate",
  seurat_clusters %in% c(2) ~ "gamma",
))

DimPlot(islet128, group.by = "manual_anno", 
        reduction = "umap", label = TRUE, repel = TRUE) +
  NoLegend()

DimPlot(islet128, group.by = "seurat_clusters", 
        reduction = "umap", label = TRUE, repel = TRUE) +
  NoLegend()


# small dotplots ----------------------------------------------------------
alpha1 <- c("GCG", "TTR", "PCSK2")
ductal1 <- c("MMP7", "KRT19", "KRT7")
gamma1 <- c("PPY", "MEIS2", "ETV1")
activated_stellate1 <- c("COL1A1", "COL1A2", "COL6A3")
beta1 <- c("INS", "IAPP", "NPTX2")
delta1 <- c("SST", "RBP4", "SEC11C")
acinar1 <- c("REG1A", "PRSS1", "PRSS2")
endothelial1 <- c("PLVAP", "RGCC", "ENG")
quiescent_stellate1 <- c("C11orf96", "CSRP2", "RGS5")
immune1 <- c("ACP5", "APOE", "C1QB")

DotPlot(islet128, features = list("beta"=beta1, "alpha" = alpha1, "delta" = delta1,
                                 "acinar" = acinar1, "ductal" = ductal1,
                                 "activated_stellate" = activated_stellate1,
                                 "endothelial" = endothelial1, "immune" = immune1, 
                                 "quiescent_stellate" = quiescent_stellate1,
                                 "gamma" = gamma1
                                 
))+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 4, 
                                   angle = 90, vjust = 0.5, hjust = 0.5))




# Heatmap -----------------------------------------------------------------

DoHeatmap(islet128, features = c(beta, alpha, delta, gamma, epsilon, cycling, ductal, 
                                endothelial, immune, quiescent_stellate,
                                schwann, activated_stellate, acinar), size = 2) +
  theme(text = element_text(size = 6))

# Variable Feature Plot ---------------------------------------------------

FeaturePlot(islet128, features = c("INS", "GCG", "SST", "PPY"), pt.size = 1)



