# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group

# first try

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))
source(here::here("R/cell_types.R"))

set.seed(100)

# islet28 load (motakis) --------------------------------------------------

# using motakis Islet28
# folder for Islet 28: /data_raw/motakis/Islet28/Solo.out/GeneFull/raw

# loading data
islet28.mtx <- ReadMtx("data_raw/motakis/Islet28/Solo.out/Gene/raw/matrix.mtx", 
                   "data_raw/motakis/Islet28/Solo.out/Gene/raw/barcodes.tsv",
                   "data_raw/motakis/Islet28/Solo.out/Gene/raw/features.tsv")

# creating plot and threshold
threshold <- valiDrops::rank_barcodes(islet28.mtx, type = "UMI") # tager virkelig lang tid!

# finding cells that pass threshold
rank.pass <- BiocGenerics::rownames(threshold$ranks[ threshold$ranks$counts >= threshold$lower.threshold, ])

# subsetting matrix
counts.subset <- islet28.mtx[, colnames(islet28.mtx) %in% rank.pass]

# creating seurat object

islet28 <- CreateSeuratObject(counts.subset, project = "islet28") #HUSK

saveRDS(islet28, file = here::here("data/seurat_objects/motakis/islet28.rds"))
islet28_raw <- readRDS(here::here("data/seurat_objects/motakis/islet28.rds"))
islet28 <- readRDS(here::here("data/seurat_objects/motakis/islet28/islet28_QC.rds"))

# Islet28 QC (motakis) -----------------------------------------------------

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

## histograms ----

### nCount_RNA
ggplot(data = islet28_raw@meta.data, aes(x=nCount_RNA)) +
  geom_histogram(bins = 200) + 
  geom_vline(xintercept = 1200)


### nFeatures_RNA
ggplot(data = islet28_raw@meta.data, aes(x=nFeature_RNA)) +
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = 800)
  

### percent.mt
ggplot(data = seurat_obj@meta.data, aes(x=percent.mt)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 15)

### exon intron ratio
ggplot(data = seurat_obj@meta.data, aes(x=exon_exonintron)) +
  geom_histogram(bins = 100)


## subsetting data ----
# remember to asses histograms when choosing thresholds
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= 800 & 
                    percent.mt < 15 &
                    nCount_RNA >= 10000 & nCount_RNA <= 40000 & 
                    exon_exonintron <= 1)

## kompleksitet af celler (nFeature/nCount) ----
islet28@meta.data$log10complexity <- log10(islet28@meta.data$nFeature_RNA)/
  log10(islet28@meta.data$nCount_RNA)

ggplot(data = islet28@meta.data, aes(x=log10complexity)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 0.80, col = "red")

## vionlinplots overview ----
VlnPlot(islet28, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## exon-intron count ----

# loader data ind som matrix (rækker = gener, kolonner = celler)
mtx_gene <- Seurat::ReadMtx(mtx = "data_raw/motakis/Islet28/Solo.out/Gene/raw/matrix.mtx",
                            cells = "data_raw/motakis/Islet28/Solo.out/Gene/raw/barcodes.tsv",
                            features = "data_raw/motakis/Islet28/Solo.out/Gene/raw/features.tsv")

mtx_genefull <- Seurat::ReadMtx(mtx = "data_raw/motakis/Islet28/Solo.out/GeneFull/raw/matrix.mtx",
                                cells = "data_raw/motakis/Islet28/Solo.out/GeneFull/raw/barcodes.tsv",
                                features = "data_raw/motakis/Islet28/Solo.out/GeneFull/raw/features.tsv")

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
seurat_obj <- Seurat::CreateSeuratObject(counts = mtx_gene_sub,
                                         assay = "RNA")

# tjekker at der er de samme celler i contrast dataframe og i seurat objektet
all.equal(rownames(contrast), rownames(seurat_obj@meta.data))

# tilføjer ratioen til meta data
seurat_obj@meta.data$exon_exonintron <- contrast$exon_exonintron

# overføre til islet28 eget data og seurat
all.equal(rownames(contrast), rownames(islet28_raw@meta.data))
# den er ens med den rå meta data - skal dette så gøres inden der subsettes??
islet28_raw@meta.data$exon_exonintron <- contrast$exon_exonintron


## normalizing data ----
islet28 <- NormalizeData(islet28)

## variable features ----
islet28 <- FindVariableFeatures(islet28, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(islet28), 10)

# plot variable features with and without labels
varfeat_plot_28 <- VariableFeaturePlot(islet28)
varfeat_plot_28_2 <- LabelPoints(plot = varfeat_plot_28, points = top10, repel = TRUE)

## scaling data and PCA ----
all.genes <- rownames(islet28)
islet28 <- ScaleData(islet28, features = all.genes)

# PCA
islet28 <- RunPCA(islet28, features = VariableFeatures(object = islet28))

# assesing important PCs
ElbowPlot(islet28) # 1:16

## marker genes ----
islet28 <- FindNeighbors(islet28, dims = 1:16)
islet28 <- FindClusters(islet28, resolution = 0.4)

head(Idents(islet28), 5)

## clustering of cells ----
islet28 <- RunUMAP(islet28, dims = 1:16)

DimPlot(islet28, reduction = "umap", label = TRUE)

# find markers
islet28.markers <- FindAllMarkers(islet28, only.pos = TRUE)

saveRDS(islet28.markers, file = here::here("data/rds_files/islet28/islet28.markers.rds"))
saveRDS(islet28, file = here::here("data/seurat_objects/motakis/islet28/islet28_QC.rds"))

islet28.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

islet28.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

FeaturePlot(islet28, features = c("ADCYAP1"),  min.cutoff = "q10")


### cluster 0 ----
cluster0.markers <- FindMarkers(islet28, ident.1 = 0, only.pos = TRUE)
head(cluster0.markers, n = 5) 
# possible beta cells (see insulin, HADH, NPTX2)
FeaturePlot(islet28, features = c("SAMD11","ADCYAP1", "MEG3", 
                                  "NPTX2", "HADH", "INS"),  min.cutoff = "q20")

# er det et enkelt cluster eller burde resolution være lavere?
cluster0.markers_comp <- FindMarkers(islet28, ident.1 = 0, ident.2 = c(10, 5), 
                                     only.pos = TRUE)
head(cluster0.markers_comp, n = 5) 
FeaturePlot(islet28, features = c("RPL34","RPS27A", "RPL24", 
                                  "RPL35A", "RPS4X"),  min.cutoff = "q20")
# generne der adskiller cluster 0 fra 10 og 5 er ikke marker gener hos Azimuth
# samt de er udtrykt hos mange clustre
# måske der burde være lavere resolution således at både 0, 5 og 10 er beta

### cluster 1 ----
cluster1.markers <- FindMarkers(islet28, ident.1 = 1)
head(cluster1.markers, n = 5) 
FeaturePlot(islet28, features = c("C5orf38","AC099509.1", "PPP1R1A", 
                                  "C12orf75", "CRYBA2"),  min.cutoff = "q20")
# måske aplha (PPP1R1A og CRYBA2)

# er det et enkelt cluster eller burde resolution være lavere?
cluster1.markers_comp <- FindMarkers(islet28, ident.1 = 1, ident.2 = c(4, 2), 
                                     only.pos = TRUE)
head(cluster1.markers_comp, n = 5) 
FeaturePlot(islet28, features = c("SERF2","TUBA1B", "COX8A", 
                                  "COX17", "C12orf75"),  min.cutoff = "q20")

# igen - resolution bør muligvis være lavere da de to andre clustre ikke skiller
# sig synderligt meget ud

### cluster 3 ----
cluster3.markers <- FindMarkers(islet28, ident.1 = 3)
head(cluster3.markers, n = 5) 
FeaturePlot(islet28, features = c("LY6H","BCHE", "CALB1", 
                                  "GPC5-AS1", "HHEX"),  min.cutoff = "q20")
# delta (HHEX og LY6H) gamma (GPC5-AS1) - skal den så rigtigt være delt op i 2?

### cluster 6 ----
cluster6.markers <- FindMarkers(islet28, ident.1 = 6)
head(cluster6.markers, n = 5) 
FeaturePlot(islet28, features = c("DEFB1","TACSTD2", "PDZK1IP1", 
                                  "KRT7", "SERPING1"),  min.cutoff = "q20")

# epsilon (DEFB1) ductal (KRT7 og SERPING1)
# nogle af dem udtrykkes også i især cluster 9 (og 8/7)

### cluster 7 ----
cluster7.markers <- FindMarkers(islet28, ident.1 = 7)
head(cluster7.markers, n = 5) 
FeaturePlot(islet28, features = c("COL6A3","COL15A1", "SFRP2", 
                                  "COL6A2", "COL3A1"),  min.cutoff = "q20")
# activated stelat (COL6A3, SFRP2 og COL3A1)
# ser rimelig isoleret ud til cluster 7 og altså ikke cluster 8 med

### cluster 8 ----
cluster8.markers <- FindMarkers(islet28, ident.1 = 8)
head(cluster8.markers, n = 5) 
FeaturePlot(islet28, features = c("ESAM","FABP4", "GPR4", 
                                  "ECSCR", "CD93"),  min.cutoff = "q20")
# quiescent_stellate (ESAM, FABP4)

### cluster 9 ----
cluster9.markers <- FindMarkers(islet28, ident.1 = 9)
head(cluster9.markers, n = 5) 
FeaturePlot(islet28, features = c("CXCL17","GSTA2", "ALDOB", 
                                  "ANPEP", "CPA2"),  min.cutoff = "q20")
# ingen matches hos Azimuth - kan læse mig til at CPA2 har noget med pancreas at gøre

### cluster 11 ----
cluster11.markers <- FindMarkers(islet28, ident.1 = 11)
head(cluster11.markers, n = 5) 
FeaturePlot(islet28, features = c("LAPTM5","TYROBP", "FCER1G", 
                                  "C1QC", "MS4A7"),  min.cutoff = "q20")
# immune (LAPTM5, TYROBP, FCER1G og C1QC)

### cluster 12 ----
cluster12.markers <- FindMarkers(islet28, ident.1 = 12)
head(cluster12.markers, n = 5) 
FeaturePlot(islet28, features = c("UBE2C","PBK", "NUSAP1", 
                                  "HJURP", "AURKB"),  min.cutoff = "q20")
# cycling (UBE2C og PBK)

### konklusion ----
# ønsker en resolution hvor cluster 0, 5 og 10 er en
# 4, 1 og 2 samles
# altså 9 clusters i alt
# ved forsøg med at sænke resolution forsvinder cluster 12 som er cyclin først
# og cluster 8/7 bliver til en 
# så jeg kan ikke få den opdeling som jeg regnede med

### cluster ID ----
# cluster 10, 0 og 5 = beta
# cluster 1, 4 og 2 = alpha
# cluster 3 = delta/gamma
# cluster 6 = epsilon/ductal
# cluster 7 = activated stelat
# cluster 8 = quiescent_stellate
# cluster 9 = ?
# cluster 11 = immune
# cluster 12 = cycling

new.cluster.ids <- c("beta_1", "alpha_1", "alpha_2", "delta/gamma", "alpha_3", 
                     "beta_2", "epsilon/ductal", "activated_stellate",
                     "quiescent_stellate", "?", "beta_3", "immune", "cycling")

names(new.cluster.ids) <- levels(islet28)
islet28 <- RenameIdents(islet28, new.cluster.ids)
DimPlot(islet28, reduction = "umap", label = TRUE, 
        label.size = 3, pt.size = 0.3, repel = TRUE) + NoLegend()


# random sample selection -------------------------------------------------

meta <- read.csv(here::here("data_raw/motakis/motakis_meta.csv"), sep = ";")
View(meta)

table(meta$disease)

samples <- meta %>%
  group_by(disease, sex) %>% # grupperer data
  slice_sample(n = 2)  %>% # udvælger tilfældigt 3 samples
  ungroup() 

print(samples %>% 
        select(donor, disease, sex))

# lav QC for alle prøverne individuelt (i hver sit script)
# gem alle plots (QC, heatmaps, dotplots, threshold, PCA (elbowplots))


# dotplots med isabell ----------------------------------------------------

DotPlot(islet28, features = list("beta"=beta, "alpha" = alpha))+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")
# lav plot med tre gener per celletype selv (lav liste med 3 gener per celletype, så dotplottet bliver delt pænt op)

DoHeatmap(islet28, features = c(beta, alpha, delta, gamma, epsilon, cycling, ductal, 
                                endothelial, immune, quiescent_stellate,
                                schwann, activated_stellate, acinar), size = 2) +
  theme(text = element_text(size = 6))
# brug heatmap som ekstra plot efter dotplottet
# nemmere at adskille og ud fra avarage

DotPlot(islet28, features = list("beta"=beta, "alpha" = alpha), 
        idents = new.cluster.ids)+
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", 
                                  high = "#A83708FF")

# senere kigger vi på en anden måde at finde markergener på

# Azimuth reference mapping -----------------------------------------------
remotes::install_github("satijalab/azimuth")
library(Azimuth)

        