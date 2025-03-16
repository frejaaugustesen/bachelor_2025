# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group
# QC of islet 44, T2Diabetic, male

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))
source(here::here("R/cell_types.R"))

set.seed(100)

islet44_raw <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet44_seurat.rds"))
islet44 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet44_seurat_QC.rds"))


# save QC  ----------------------------------------------------------------
# function to save QC as i go
saveRDS(islet44, file = here::here("data/seurat_objects/motakis/selected_samples/islet44_seurat_QC.rds"))


# load data ---------------------------------------------------------------

# loader data ind som matrix (rækker = gener, kolonner = celler)
mtx_gene <- Seurat::ReadMtx(mtx = "data_raw/motakis/Islet44/Solo.out/Gene/raw/matrix.mtx",
                            cells = "data_raw/motakis/Islet44/Solo.out/Gene/raw/barcodes.tsv",
                            features = "data_raw/motakis/Islet44/Solo.out/Gene/raw/features.tsv")

mtx_genefull <- Seurat::ReadMtx(mtx = "data_raw/motakis/Islet44/Solo.out/GeneFull/raw/matrix.mtx",
                                cells = "data_raw/motakis/Islet44/Solo.out/GeneFull/raw/barcodes.tsv",
                                features = "data_raw/motakis/Islet44/Solo.out/GeneFull/raw/features.tsv")

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
islet44 <- Seurat::CreateSeuratObject(counts = mtx_gene_sub,
                                      assay = "RNA",
                                      project = "islet44")

# tjekker at der er de samme celler i contrast dataframe og i seurat objektet
all.equal(rownames(contrast), rownames(islet44@meta.data))

# tilføjer ratioen til meta data
islet44@meta.data$exon_exonintron <- contrast$exon_exonintron

# gemmer seuratobjekt
saveRDS(islet44, file = here::here("data/seurat_objects/motakis/selected_samples/islet44_seurat.rds"))


# QC ----------------------------------------------------------------------

# adding percentage of mitochondrial RNA ----
islet44[["percent.mt"]] <- PercentageFeatureSet(islet44, pattern = "^MT-")

# kompleksitet af celler (nFeature/nCount) ----
islet44@meta.data$log10complexity <- log10(islet44@meta.data$nFeature_RNA)/
  log10(islet44@meta.data$nCount_RNA)

## histograms ----

### nCount_RNA
ggplot(data = islet44@meta.data, aes(x=nCount_RNA)) +
  geom_histogram(bins = 200) + 
  geom_vline(xintercept = 7000) +
  geom_vline(xintercept = 80000)


### nFeatures_RNA
ggplot(data = islet44@meta.data, aes(x=nFeature_RNA)) +
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = 2000)


### percent.mt
ggplot(data = islet44@meta.data, aes(x=percent.mt)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 15)

### exon intron ratio
ggplot(data = islet44@meta.data, aes(x=exon_exonintron)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 1)

### complexity
ggplot(data = islet44@meta.data, aes(x=log10complexity)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 0.80)

## subsetting data ----
# remember to asses histograms when choosing thresholds
islet44 <- subset(islet44, subset = nFeature_RNA >= 2000 & 
                    percent.mt < 15 &
                    nCount_RNA >= 7000 & nCount_RNA <= 80000 & 
                    exon_exonintron <= 1 & log10complexity >= 0.80)

# saving QC as seurat object - to get the previous histograms, load raw seurat
saveRDS(islet44, file = here::here("data/seurat_objects/motakis/selected_samples/islet44_seurat_QC.rds"))

## normalizing data ----
islet44 <- NormalizeData(islet44)

## variable features ----
islet44 <- FindVariableFeatures(islet44,
                                selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(islet44), 10)

# plot variable features with and without labels
varfeat_plot_44 <- VariableFeaturePlot(islet44)
varfeat_plot_44_2 <- LabelPoints(plot = varfeat_plot_44,
                                 points = top10, repel = TRUE)

## scaling data and PCA ----
all.genes <- rownames(islet44)
islet44 <- ScaleData(islet44, features = all.genes)

# PCA
islet44 <- RunPCA(islet44, features = VariableFeatures(object = islet44))

# assesing important PCs
ElbowPlot(islet44) # 1:10

# clustering  -------------------------------------------------------------

## marker genes ----
islet44 <- FindNeighbors(islet44, dims = 1:10)
islet44 <- FindClusters(islet44, resolution = 0.6)

## clustering of cells ----
islet44 <- RunUMAP(islet44, dims = 1:10)

DimPlot(islet44, reduction = "umap", label = TRUE)

# Dotplots ----------------------------------------------------------------

DotPlot1 <- DotPlot(islet44, 
                    features = list("beta"=beta, "alpha" = alpha, 
                                    "delta" = delta, "gamma" = gamma
                    ))+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3.5))

DotPlot2 <- DotPlot(islet44, 
                    features = list("acinar" = acinar, "ductal" = ductal, 
                                    "cycling" = cycling, "immune" = immune
                    ))+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3.5))

DotPlot3 <- DotPlot(islet44, 
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
# cluster 1: alpha
# cluster 2: delta
# cluster 3: beta/gamma
# cluster 4: ductal
# cluster 5: ductal
# cluster 6: activated stellate
# cluster 7: alpha
# cluster 8: acinar
# cluster 9: quiescent stellate/endothelial
# cluster 10: immune


## klade til annotering ----

islet44@meta.data <- islet44@meta.data %>% dplyr::mutate(manual_anno = dplyr::case_when(
  seurat_clusters %in% c(0, 1, 7) ~ "alpha",
  seurat_clusters %in% c(3) ~ "beta/gamma",
  seurat_clusters %in% c(6) ~ "activated_stellate",
  seurat_clusters %in% c(8) ~ "acinar",
  seurat_clusters %in% c(4, 5) ~ "ductal",
  seurat_clusters %in% c(2) ~ "delta",
  seurat_clusters %in% c(9) ~ "endothelial/quiescent_stellate",
  seurat_clusters %in% c(10) ~ "immune",
))

DimPlot(islet44, group.by = "manual_anno", 
        reduction = "umap", label = TRUE, repel = TRUE) +
  NoLegend()

DimPlot(islet44, group.by = "seurat_clusters", 
        reduction = "umap", label = TRUE, repel = TRUE) +
  NoLegend()


# small dotplots ----------------------------------------------------------
alpha1 <- c("GCG", "TTR", "PCSK2")
beta1 <- c("INS", "IAPP", "NPTX2")
delta1 <- c("SST", "RBP4", "SEC11C")
gamma1 <- c("PPY", "MEIS2", "ETV1")
ductal1 <- c("ANXA4", "KRT19", "KRT7")
activated_stellate1 <- c("COL1A1", "COL1A2", "COL6A3")
acinar1 <- c("REG1A", "PRSS1", "PRSS2")
endothelial1 <- c("PLVAP", "RGCC", "ENG")
immune1 <- c("ACP5", "APOE", "C1QB")
quiescent_stellate1 <- c("C11orf96", "CSRP2", "ESAM")

DotPlot(islet44, features = list("beta"=beta1, "alpha" = alpha1, "delta" = delta1,
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

DoHeatmap(islet44, features = c(beta, alpha, delta, gamma, epsilon, cycling, ductal, 
                                endothelial, immune, quiescent_stellate,
                                schwann, activated_stellate, acinar), size = 2) +
  theme(text = element_text(size = 6))


# Variable Feature Plot ---------------------------------------------------

FeaturePlot(islet44, features = c("INS", "GCG", "SST", "PPY"), pt.size = 1)



