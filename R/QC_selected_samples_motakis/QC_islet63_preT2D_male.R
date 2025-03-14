# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group
# QC of islet 63, Pre-Diabetic, male

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))
source(here::here("R/cell_types.R"))

set.seed(100)

islet63_raw <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet63_seurat.rds"))
islet63 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet63_seurat_QC.rds"))


# save QC  ----------------------------------------------------------------
# function to save QC as i go
saveRDS(islet63, file = here::here("data/seurat_objects/motakis/selected_samples/islet63_seurat_QC.rds"))


# load data ---------------------------------------------------------------

# loader data ind som matrix (rækker = gener, kolonner = celler)
mtx_gene <- Seurat::ReadMtx(mtx = "data_raw/motakis/Islet63/Solo.out/Gene/raw/matrix.mtx",
                            cells = "data_raw/motakis/Islet63/Solo.out/Gene/raw/barcodes.tsv",
                            features = "data_raw/motakis/Islet63/Solo.out/Gene/raw/features.tsv")

mtx_genefull <- Seurat::ReadMtx(mtx = "data_raw/motakis/Islet63/Solo.out/GeneFull/raw/matrix.mtx",
                                cells = "data_raw/motakis/Islet63/Solo.out/GeneFull/raw/barcodes.tsv",
                                features = "data_raw/motakis/Islet63/Solo.out/GeneFull/raw/features.tsv")

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
islet63 <- Seurat::CreateSeuratObject(counts = mtx_gene_sub,
                                      assay = "RNA",
                                      project = "islet63")

# tjekker at der er de samme celler i contrast dataframe og i seurat objektet
all.equal(rownames(contrast), rownames(islet63@meta.data))

# tilføjer ratioen til meta data
islet63@meta.data$exon_exonintron <- contrast$exon_exonintron

# gemmer seuratobjekt
saveRDS(islet63, file = here::here("data/seurat_objects/motakis/selected_samples/islet63_seurat.rds"))


# QC ----------------------------------------------------------------------

# adding percentage of mitochondrial RNA ----
islet63[["percent.mt"]] <- PercentageFeatureSet(islet63, pattern = "^MT-")

# kompleksitet af celler (nFeature/nCount) ----
islet63@meta.data$log10complexity <- log10(islet63@meta.data$nFeature_RNA)/
  log10(islet63@meta.data$nCount_RNA)

## histograms ----

### nCount_RNA
ggplot(data = islet63@meta.data, aes(x=nCount_RNA)) +
  geom_histogram(bins = 200) + 
  geom_vline(xintercept = 11800) +
  geom_vline(xintercept = 40000)


### nFeatures_RNA
ggplot(data = islet63@meta.data, aes(x=nFeature_RNA)) +
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = 2000)


### percent.mt
ggplot(data = islet63@meta.data, aes(x=percent.mt)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 20)

### exon intron ratio
ggplot(data = islet63@meta.data, aes(x=exon_exonintron)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 1)

### complexity
ggplot(data = islet63@meta.data, aes(x=log10complexity)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 0.80)

## subsetting data ----
# remember to asses histograms when choosing thresholds
islet63 <- subset(islet63, subset = nFeature_RNA >= 2000 & 
                    percent.mt < 20 &
                    nCount_RNA >= 11800 & nCount_RNA <= 40000 & 
                    exon_exonintron <= 1 & log10complexity >= 0.80)

# saving QC as seurat object - to get the previous histograms, load raw seurat
saveRDS(islet63, file = here::here("data/seurat_objects/motakis/selected_samples/islet63_seurat_QC.rds"))

## normalizing data ----
islet63 <- NormalizeData(islet63)

## variable features ----
islet63 <- FindVariableFeatures(islet63,
                                selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(islet63), 10)

# plot variable features with and without labels
varfeat_plot_63 <- VariableFeaturePlot(islet63)
varfeat_plot_63_2 <- LabelPoints(plot = varfeat_plot_63,
                                 points = top10, repel = TRUE)

## scaling data and PCA ----
all.genes <- rownames(islet63)
islet63 <- ScaleData(islet63, features = all.genes)

# PCA
islet63 <- RunPCA(islet63, features = VariableFeatures(object = islet63))

# assesing important PCs
ElbowPlot(islet63) # 1:16

# clustering  -------------------------------------------------------------

## marker genes ----
islet63 <- FindNeighbors(islet63, dims = 1:16)
islet63 <- FindClusters(islet63, resolution = 0.6)

## clustering of cells ----
islet63 <- RunUMAP(islet63, dims = 1:16)

DimPlot(islet63, reduction = "umap", label = TRUE)

# Dotplots ----------------------------------------------------------------

DotPlot1 <- DotPlot(islet63, 
                    features = list("beta"=beta, "alpha" = alpha, 
                                    "delta" = delta, "gamma" = gamma
                    ))+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3.5))

DotPlot2 <- DotPlot(islet63, 
                    features = list("acinar" = acinar, "ductal" = ductal, 
                                    "cycling" = cycling, "immune" = immune
                    ))+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3.5))

DotPlot3 <- DotPlot(islet63, 
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

# cluster 0: acinar
# cluster 1: alpha
# cluster 2: ductal
# cluster 3: acinar
# cluster 4: delta
# cluster 5: activated stellate
# cluster 6: beta
# cluster 7: immune


## klade til annotering ----

new.cluster.ids <- c("acinar", "alpha", "ductal", "acinar", "delta",
                     "activated_stellate", "beta", "immune")


names(new.cluster.ids) <- levels(islet63)
islet63 <- RenameIdents(islet63, new.cluster.ids)
DimPlot(islet63, reduction = "umap", label = TRUE, 
        label.size = 3, pt.size = 0.5, repel = TRUE) + NoLegend()


# small dotplots ----------------------------------------------------------
acinar1 <- c("REG1B", "CELA3A", "CPA1")
alpha1 <- c("GCG", "GC", "PCSK2")
ductal1 <- c("MMP7", "CFTR", "ANXA4")
delta1 <- c("SST", "RBP4", "LEPR")
activated_stellate1 <- c("COL1A2", "COL3A1", "BGN")
beta1 <- c("IAPP", "INS", "NPTX2")
immune1 <- c("ACP5", "APOE", "C1QB")

DotPlot(islet63, features = list("beta"=beta1, "alpha" = alpha1, "delta" = delta1,
                                 "acinar" = acinar1, "ductal" = ductal1,
                                 "activated_stellate" = activated_stellate1,
                                 "immune" = immune1 
                                 
))+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 4, 
                                   angle = 90, vjust = 0.5, hjust = 0.5))




# Heatmap -----------------------------------------------------------------

DoHeatmap(islet63, features = c(beta, alpha, delta, gamma, epsilon, cycling, ductal, 
                                endothelial, immune, quiescent_stellate,
                                schwann, activated_stellate, acinar), size = 2) +
  theme(text = element_text(size = 6))


