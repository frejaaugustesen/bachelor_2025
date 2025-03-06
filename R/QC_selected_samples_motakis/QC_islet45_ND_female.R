# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group
# QC of islet 45, NonDiabetic, female

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))
source(here::here("R/cell_types.R"))

set.seed(100)

islet45_raw <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet45_seurat.rds"))
islet45 <- readRDS(here::here("data/seurat_objects/motakis/selected_samples/islet45_seurat_QC.rds"))


# save QC  ----------------------------------------------------------------
# function to save QC as i go
saveRDS(islet45, file = here::here("data/seurat_objects/motakis/selected_samples/islet45_seurat_QC.rds"))




# load data ---------------------------------------------------------------

# loader data ind som matrix (rækker = gener, kolonner = celler)
mtx_gene <- Seurat::ReadMtx(mtx = "data_raw/motakis/Islet45/Solo.out/Gene/raw/matrix.mtx",
                            cells = "data_raw/motakis/Islet45/Solo.out/Gene/raw/barcodes.tsv",
                            features = "data_raw/motakis/Islet45/Solo.out/Gene/raw/features.tsv")

mtx_genefull <- Seurat::ReadMtx(mtx = "data_raw/motakis/Islet45/Solo.out/GeneFull/raw/matrix.mtx",
                                cells = "data_raw/motakis/Islet45/Solo.out/GeneFull/raw/barcodes.tsv",
                                features = "data_raw/motakis/Islet45/Solo.out/GeneFull/raw/features.tsv")

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
islet45 <- Seurat::CreateSeuratObject(counts = mtx_gene_sub,
                                      assay = "RNA")

# tjekker at der er de samme celler i contrast dataframe og i seurat objektet
all.equal(rownames(contrast), rownames(islet45@meta.data))

# tilføjer ratioen til meta data
islet45@meta.data$exon_exonintron <- contrast$exon_exonintron

# gemmer seuratobjekt
saveRDS(islet45, file = here::here("data/seurat_objects/motakis/selected_samples/islet45_seurat.rds"))


# QC ----------------------------------------------------------------------

# adding percentage of mitochondrial RNA ----
islet45[["percent.mt"]] <- PercentageFeatureSet(islet45, pattern = "^MT-")

# kompleksitet af celler (nFeature/nCount) ----
islet45@meta.data$log10complexity <- log10(islet45@meta.data$nFeature_RNA)/
  log10(islet45@meta.data$nCount_RNA)

## histograms ----

### nCount_RNA
ggplot(data = islet45@meta.data, aes(x=nCount_RNA)) +
  geom_histogram(bins = 200) + 
  geom_vline(xintercept = 4000) +
  geom_vline(xintercept = 40000)


### nFeatures_RNA
ggplot(data = islet45@meta.data, aes(x=nFeature_RNA)) +
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = 1200)


### percent.mt
ggplot(data = islet45@meta.data, aes(x=percent.mt)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 25)

### exon intron ratio
ggplot(data = islet45@meta.data, aes(x=exon_exonintron)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 1)

### complexity
ggplot(data = islet45@meta.data, aes(x=log10complexity)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 0.80)

## subsetting data ----
# remember to asses histograms when choosing thresholds
islet45 <- subset(islet45, subset = nFeature_RNA >= 1200 & 
                    percent.mt < 25 &
                    nCount_RNA >= 4000 & nCount_RNA <= 40000 & 
                    exon_exonintron <= 1)

# saving QC as seurat object - to get the previous histograms, load raw seurat
saveRDS(islet45, file = here::here("data/seurat_objects/motakis/selected_samples/islet45_seurat_QC.rds"))

## normalizing data ----
islet45 <- NormalizeData(islet45)

## variable features ----
islet45 <- FindVariableFeatures(islet45,
                                selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(islet45), 10)

# plot variable features with and without labels
varfeat_plot_45 <- VariableFeaturePlot(islet45)
varfeat_plot_45_2 <- LabelPoints(plot = varfeat_plot_45,
                                 points = top10, repel = TRUE)

## scaling data and PCA ----
all.genes <- rownames(islet45)
islet45 <- ScaleData(islet45, features = all.genes)

# PCA
islet45 <- RunPCA(islet45, features = VariableFeatures(object = islet45))

# assesing important PCs
ElbowPlot(islet45) # 1:15

