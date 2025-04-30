# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))

set.seed(100)

# load --------------------------------------------------------------------

wang <- qread("/work/bachelor_2025/data/seurat_objects/wang_seurat_QC.qs")
motakis <- qread("/work/bachelor_2025/data/seurat_objects/motakis/islet_all_subset.qs")

wang_beta <- qread("/work/bachelor_2025/data/seurat_objects/wang_beta.qs")
motakis_beta <- qread("/work/bachelor_2025/data/seurat_objects/motakis_beta.qs")

wang_beta_sub <- qread("/work/bachelor_2025/data/seurat_objects/wang_beta_sub.qs")
wang_beta_train <- qread("/work/bachelor_2025/data/seurat_objects/wang_beta_train.qs")
wang_beta_test <- qread("/work/bachelor_2025/data/seurat_objects/wang_beta_test.qs")

var_genes <- qread("/work/bachelor_2025/data/var_genes.qs")

bst <- qread(here::here("data/xgboost/finished_model/bst.final.qs"))
wrong <- qread(here::here("data/xgboost/finished_model/wrong_vector.qs"))

# saving objects ----
qsave(motakis_beta, file = here::here("data/seurat_objects/motakis_beta_integrated.qs"))

# correcting motakis data --------------------------------------------------
# need to scale and normalize

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

motakis_beta <- FindNeighbors(motakis_beta, dims = 1:15, reduction = "pca")
motakis_beta <- FindClusters(motakis_beta)

motakis_beta <- RunUMAP(motakis_beta, dims = 1:15)

DimPlot(motakis_beta, reduction = "umap", label = TRUE)

DimPlot(motakis_beta, reduction = "umap", label = TRUE, repel = TRUE,
        label.size = 3,
        group.by = "manual_anno") +
  NoLegend()

# integration
motakis_beta <- IntegrateLayers(object = motakis_beta, 
                             method = HarmonyIntegration, 
                             orig.reduction = "pca", 
                             new.reduction = "HarmonyIntegration",
                             verbose = FALSE)

# re-join layers after integration
motakis_beta[["RNA"]] <- JoinLayers(motakis_beta[["RNA"]])

motakis_beta <- FindNeighbors(motakis_beta, reduction = "HarmonyIntegration", dims = 1:15)
motakis_beta <- FindClusters(motakis_beta)

motakis_beta <- RunUMAP(motakis_beta, dims = 1:15, reduction = "HarmonyIntegration")


# correct attempt ----------------------------------------------------------
# i will try to attempt to recreate code from the wang sander study
# https://github.com/gaoweiwang/Islet_snATACseq/blob/8145d545d979507193fd5bef71b1b4119bed6e66/scripts/AI_subtype.ipynb#L710

## setting up data ----

# creating data_use
# i already have the matrix for train and test in such a matrix...
# this should only be one matrix as i am doing leave one out
# using only var_genes
# Get the log-normalized data matrix - not using the sub one as we sort that out later
rna_data <- GetAssayData(motakis_beta, assay = "RNA", slot = "data")

# keeping genes, which are in var_genes (determined earlier)
data_use <- rna_data[var_genes, ]

# assesing data_use
length(data_use)
rownames(data_use)
colnames(data_use)

# changing values to only being either 0 or 1
data_use@x[data_use@x > 0] <- 1

# now i have a matrix with the rownames being the genes and the 
# colnames being donorid_barcode with either 1 or 0
# this needs to be changed so rows becomes columns and columns becomes rows
data_use <- t(data_use)

# creating vector with all donor names for ND and T2D
motakis_beta_sub <- subset(motakis_beta, subset = disease != "pre") # used earlier
motakis_meta_T2D_ND <- motakis_beta_sub@meta.data
donor_all <- unique(motakis_meta_T2D_ND$orig.ident)

## training model (loop) ------------------------------------------------------

# first creating dataframe with disease and predicted disease 
M_new <- motakis_beta@meta.data #husk at tilføje disease længere oppe
M_new$pre_disease <- motakis_beta@meta.data$disease
M_new$subtype <- motakis_beta@meta.data$disease
M_new <- dplyr::rename(M_new, donor = orig.ident)

# renaming disease states
M_new <- M_new %>%
  mutate(subtype = recode(subtype,
                          "NonDiabetic" = "nd",
                          "PreDiabetic" = "pre",
                          "Type2Diabetic" = "t2d"))

M_new <- M_new %>%
  mutate(disease = recode(disease,
                          "NonDiabetic" = "nd",
                          "PreDiabetic" = "pre",
                          "Type2Diabetic" = "t2d"))

M_new <- M_new %>%
  mutate(pre_disease = recode(pre_disease,
                          "NonDiabetic" = "nd",
                          "PreDiabetic" = "pre",
                          "Type2Diabetic" = "t2d"))

# Predict på prediabetes (og de donorer der trænet på) ----
keep=as.character(M_new$disease) %in% c("t2d", "nd", "pre")
pre.x = data_use[keep,]
pre.x1=as(pre.x, "dgCMatrix")
pred <- predict(bst, pre.x1)

rowSums(table(pred > 0.5, M_new[keep, "donor"]))
table(pred > 0.5, M_new[keep, "donor"])
# false er nd og true t2d

## adding prediction to meta data ------------------------------------
# nu vil jeg merge predictions fra M_new (for prediabetics) over i metadata
# alt dette er fra loopet, jeg har genbrugt
temp_label=rep('no',dim(M_new)[1])
temp_pro=rep(-1,dim(M_new)[1])

pred_label=as.numeric(pred > 0.5)
pred_label1=pred_label
pred_label1[pred_label==0]='nd' 
pred_label1[pred_label==1]='t2d' 

keep_test=as.character(M_new$disease) %in% c("t2d", "nd", "pre")

temp_label[keep_test]=pred_label1
temp_pro[keep_test]=pred

M_new_pre <- M_new # denne er ny, for ikke at ødelægge originalen

M_new_pre$subtype=temp_label 
M_new_pre$subtype_probability=temp_pro
View(M_new_pre)

M_new_pre <-  M_new_pre %>% 
  tibble::rownames_to_column("barcode")

# nu har jeg en dataframe som jeg kan bruge til mit meta data i mit pre seurat objekt
motakis_beta@meta.data <-  motakis_beta@meta.data %>% 
  tibble::rownames_to_column("barcode") %>% 
  dplyr::left_join(y = M_new_pre %>% select(barcode, subtype, subtype_probability),
                   by = "barcode") %>% 
  tibble::column_to_rownames("barcode")

motakis_beta@meta.data <- motakis_beta@meta.data %>%
  mutate(disease = recode(disease,
                          "NonDiabetic" = "nd",
                          "PreDiabetic" = "pre",
                          "Type2Diabetic" = "t2d"))

# plot of subtypes distribution -------------------------------------------

# getting an overview
table(motakis_beta@meta.data$subtype)
table(motakis_beta@meta.data$disease)
correct <- sum(motakis_beta@meta.data$subtype == motakis_beta@meta.data$disease)
total <- length(motakis_beta@meta.data$subtype)

correct_percent <- (correct/total)*100

motakis_beta@meta.data %>%
  mutate(match = subtype == disease) %>%
  group_by(disease, match) %>%
  summarise(n = n(), .groups = "drop")


ggplot(motakis_beta@meta.data, aes(x = orig.ident, fill = subtype)) +
  geom_bar(position = "fill") +
  facet_wrap(~ disease) +
  labs(
    title = "Subtype distribution in prediabetics",
    x = "Donor",
    y = "Percentage"
  ) +
  theme(axis.text.x = element_text(angle = 45))

motakis_beta@meta.data %>%
  filter(orig.ident != "" & !is.na(orig.ident)) %>%
  group_by(disease) %>%
  filter(orig.ident %in% unique(orig.ident)) %>%  # optional if filtering earlier
  mutate(orig.ident = factor(orig.ident)) %>%
  ungroup() %>%
  ggplot(aes(x = orig.ident, fill = subtype)) +
  geom_bar(position = "fill") +
  facet_wrap(~ disease, scales = "free_x") +  # this is the real key: free x-axis per facet
  labs(
    title = "Subtype distribution by donor and disease group",
    x = "Donor",
    y = "Percentage"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

motakis_beta@meta.data %>%
  group_by(disease) %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  ungroup() %>%
  ggplot(aes(x = orig.ident, fill = subtype)) +
  geom_bar(position = "fill") +
  facet_wrap(~ disease, scales = "free_x") +
  labs(
    title = "Predicted subtypes per donor",
    x = "Donor",
    y = "Percentage"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# clustering of subtypes --------------------------------------------------
DimPlot(motakis_beta, reduction = "umap", group.by = "subtype")
