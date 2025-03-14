# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group
# QC of islet 52, T2Diabetic, female

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))
source(here::here("R/cell_types.R"))

set.seed(100)

islet52_raw <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet52_seurat.rds"))
islet52 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet52_seurat_QC.rds"))


# save QC  ----------------------------------------------------------------
# function to save QC as i go
saveRDS(islet52, file = here::here("data/seurat_objects/motakis/selected_samples/islet52_seurat_QC.rds"))


# load data ---------------------------------------------------------------

# loader data ind som matrix (rækker = gener, kolonner = celler)
mtx_gene <- Seurat::ReadMtx(mtx = "data_raw/motakis/Islet52/Solo.out/Gene/raw/matrix.mtx",
                            cells = "data_raw/motakis/Islet52/Solo.out/Gene/raw/barcodes.tsv",
                            features = "data_raw/motakis/Islet52/Solo.out/Gene/raw/features.tsv")

mtx_genefull <- Seurat::ReadMtx(mtx = "data_raw/motakis/Islet52/Solo.out/GeneFull/raw/matrix.mtx",
                                cells = "data_raw/motakis/Islet52/Solo.out/GeneFull/raw/barcodes.tsv",
                                features = "data_raw/motakis/Islet52/Solo.out/GeneFull/raw/features.tsv")

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
islet52 <- Seurat::CreateSeuratObject(counts = mtx_gene_sub,
                                      assay = "RNA",
                                      project = "islet52")

# tjekker at der er de samme celler i contrast dataframe og i seurat objektet
all.equal(rownames(contrast), rownames(islet52@meta.data))

# tilføjer ratioen til meta data
islet52@meta.data$exon_exonintron <- contrast$exon_exonintron

# gemmer seuratobjekt
saveRDS(islet52, file = here::here("data/seurat_objects/motakis/selected_samples/islet52_seurat.rds"))


# QC ----------------------------------------------------------------------

# adding percentage of mitochondrial RNA ----
islet52[["percent.mt"]] <- PercentageFeatureSet(islet52, pattern = "^MT-")

# kompleksitet af celler (nFeature/nCount) ----
islet52@meta.data$log10complexity <- log10(islet52@meta.data$nFeature_RNA)/
  log10(islet52@meta.data$nCount_RNA)

## histograms ----

### nCount_RNA
ggplot(data = islet52@meta.data, aes(x=nCount_RNA)) +
  geom_histogram(bins = 200) + 
  geom_vline(xintercept = 3000) +
  geom_vline(xintercept = 30000)


### nFeatures_RNA
ggplot(data = islet52@meta.data, aes(x=nFeature_RNA)) +
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = 1200)


### percent.mt
ggplot(data = islet52@meta.data, aes(x=percent.mt)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 15)

### exon intron ratio
ggplot(data = islet52@meta.data, aes(x=exon_exonintron)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 1)

### complexity
ggplot(data = islet52@meta.data, aes(x=log10complexity)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 0.80)

## subsetting data ----
# remember to asses histograms when choosing thresholds
islet52 <- subset(islet52, subset = nFeature_RNA >= 1200 & 
                    percent.mt < 15 &
                    nCount_RNA >= 3000 & nCount_RNA <= 30000 & 
                    exon_exonintron <= 1 & log10complexity >= 0.80)

# saving QC as seurat object - to get the previous histograms, load raw seurat
saveRDS(islet52, file = here::here("data/seurat_objects/motakis/selected_samples/islet52_seurat_QC.rds"))

## normalizing data ----
islet52 <- NormalizeData(islet52)

## variable features ----
islet52 <- FindVariableFeatures(islet52,
                                selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(islet52), 10)

# plot variable features with and without labels
varfeat_plot_52 <- VariableFeaturePlot(islet52)
varfeat_plot_52_2 <- LabelPoints(plot = varfeat_plot_52,
                                 points = top10, repel = TRUE)

## scaling data and PCA ----
all.genes <- rownames(islet52)
islet52 <- ScaleData(islet52, features = all.genes)

# PCA
islet52 <- RunPCA(islet52, features = VariableFeatures(object = islet52))

# assesing important PCs
ElbowPlot(islet52) # 1:16

# clustering  -------------------------------------------------------------

## marker genes ----
islet52 <- FindNeighbors(islet52, dims = 1:16)
islet52 <- FindClusters(islet52, resolution = 0.6)

## clustering of cells ----
islet52 <- RunUMAP(islet52, dims = 1:16)

DimPlot(islet52, reduction = "umap", label = TRUE)

# Dotplots ----------------------------------------------------------------

DotPlot1 <- DotPlot(islet52, 
                    features = list("beta"=beta, "alpha" = alpha, 
                                    "delta" = delta, "gamma" = gamma
                    ))+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3.5))

DotPlot2 <- DotPlot(islet52, 
                    features = list("acinar" = acinar, "ductal" = ductal, 
                                    "cycling" = cycling, "immune" = immune
                    ))+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3.5))

DotPlot3 <- DotPlot(islet52, 
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
# cluster 2: activated stellate
# cluster 3: ductal
# cluster 4: beta/delta
# cluster 5: endothelial
# cluster 6: alpha
# cluster 7: acinar
# cluster 8: ductal
# cluster 9: quiescent stellate
# cluster 10: immune


## klade til annotering ----

new.cluster.ids <- c("alpha", "alpha", "activated_stellate", "ductal", 
                     "beta/delta", "endothelial", "alpha", "acinar",
                     "ductal", "quiescent_stellate", "immune")


names(new.cluster.ids) <- levels(islet52)
islet52 <- RenameIdents(islet52, new.cluster.ids)
DimPlot(islet52, reduction = "umap", label = TRUE, 
        label.size = 3, pt.size = 0.5, repel = TRUE) + NoLegend()


# small dotplots ----------------------------------------------------------
alpha1 <- c("GCG", "TTR", "PCSK2")
beta1 <- c("INS", "IAPP", "HADH")
ductal1 <- c("SPP1", "ANXA4", "KRT19")
activated_stellate1 <- c("COL1A1", "COL1A2", "COL6A3")
acinar1 <- c("REG1A", "PRSS1", "PRSS2")
delta1 <- c("SST", "RBP4", "SEC11C")
endothelial1 <- c("PLVAP", "RGCC", "ENG")
immune1 <- c("ACP5", "APOE", "C1QB")
quiescent_stellate1 <- c("C11orf96", "CSRP2", "RGS5")

DotPlot(islet52, features = list("beta"=beta1, "alpha" = alpha1, "delta" = delta1,
                                 "acinar" = acinar1, "ductal" = ductal1,
                                 "activated_stellate" = activated_stellate1,
                                "endothelial" = endothelial1, "immune" = immune1, 
                                "quiescent_stellate" = quiescent_stellate1
                                 
))+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")+
  theme(text = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 4, 
                                   angle = 90, vjust = 0.5, hjust = 0.5))



# Heatmap -----------------------------------------------------------------

DoHeatmap(islet52, features = c(beta, alpha, delta, gamma, epsilon, cycling, ductal, 
                                endothelial, immune, quiescent_stellate,
                                schwann, activated_stellate, acinar), size = 2) +
  theme(text = element_text(size = 6))
