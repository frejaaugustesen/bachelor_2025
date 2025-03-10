# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group
# QC of islet 127, Pre-Diabetic, male

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))
source(here::here("R/cell_types.R"))

set.seed(100)

islet127_raw <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet127_seurat.rds"))
islet127 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet127_seurat_QC.rds"))


# save QC  ----------------------------------------------------------------
# function to save QC as i go
saveRDS(islet127, file = here::here("data/seurat_objects/motakis/selected_samples/islet127_seurat_QC.rds"))


# load data ---------------------------------------------------------------

# loader data ind som matrix (rækker = gener, kolonner = celler)
mtx_gene <- Seurat::ReadMtx(mtx = "data_raw/motakis/Islet127/Solo.out/Gene/raw/matrix.mtx",
                            cells = "data_raw/motakis/Islet127/Solo.out/Gene/raw/barcodes.tsv",
                            features = "data_raw/motakis/Islet127/Solo.out/Gene/raw/features.tsv")

mtx_genefull <- Seurat::ReadMtx(mtx = "data_raw/motakis/Islet127/Solo.out/GeneFull/raw/matrix.mtx",
                                cells = "data_raw/motakis/Islet127/Solo.out/GeneFull/raw/barcodes.tsv",
                                features = "data_raw/motakis/Islet127/Solo.out/GeneFull/raw/features.tsv")

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
islet127 <- Seurat::CreateSeuratObject(counts = mtx_gene_sub,
                                      assay = "RNA")

# tjekker at der er de samme celler i contrast dataframe og i seurat objektet
all.equal(rownames(contrast), rownames(islet127@meta.data))

# tilføjer ratioen til meta data
islet127@meta.data$exon_exonintron <- contrast$exon_exonintron

# gemmer seuratobjekt
saveRDS(islet127, file = here::here("data/seurat_objects/motakis/selected_samples/islet127_seurat.rds"))


# QC ----------------------------------------------------------------------

# adding percentage of mitochondrial RNA ----
islet127[["percent.mt"]] <- PercentageFeatureSet(islet127, pattern = "^MT-")

# kompleksitet af celler (nFeature/nCount) ----
islet127@meta.data$log10complexity <- log10(islet127@meta.data$nFeature_RNA)/
  log10(islet127@meta.data$nCount_RNA)

## histograms ----

### nCount_RNA
ggplot(data = islet127@meta.data, aes(x=nCount_RNA)) +
  geom_histogram(bins = 200) + 
  geom_vline(xintercept = 12000) +
  geom_vline(xintercept = 60000)


### nFeatures_RNA
ggplot(data = islet127@meta.data, aes(x=nFeature_RNA)) +
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = 1700)


### percent.mt
ggplot(data = islet127@meta.data, aes(x=percent.mt)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 20)

### exon intron ratio
ggplot(data = islet127@meta.data, aes(x=exon_exonintron)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 1)

### complexity
ggplot(data = islet127@meta.data, aes(x=log10complexity)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 0.80)

## subsetting data ----
# remember to asses histograms when choosing thresholds
islet127 <- subset(islet127, subset = nFeature_RNA >= 1700 & 
                    percent.mt < 20 &
                    nCount_RNA >= 12000 & nCount_RNA <= 60000 & 
                    exon_exonintron <= 1)

# saving QC as seurat object - to get the previous histograms, load raw seurat
saveRDS(islet127, file = here::here("data/seurat_objects/motakis/selected_samples/islet127_seurat_QC.rds"))

## normalizing data ----
islet127 <- NormalizeData(islet127)

## variable features ----
islet127 <- FindVariableFeatures(islet127,
                                selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(islet127), 10)

# plot variable features with and without labels
varfeat_plot_127 <- VariableFeaturePlot(islet127)
varfeat_plot_127_2 <- LabelPoints(plot = varfeat_plot_127,
                                 points = top10, repel = TRUE)

## scaling data and PCA ----
all.genes <- rownames(islet127)
islet127 <- ScaleData(islet127, features = all.genes)

# PCA
islet127 <- RunPCA(islet127, features = VariableFeatures(object = islet127))

# assesing important PCs
ElbowPlot(islet127) # 1:16
