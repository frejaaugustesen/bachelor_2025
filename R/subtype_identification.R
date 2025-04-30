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

matrix_data_train <- qread("/work/bachelor_2025/data/matrix_data_train.qs")
matrix_data_test <- qread("/work/bachelor_2025/data/matrix_data_test.qs")

bst <- qread(here::here("data/xgboost/finished_model/bst.final.qs"))
wrong <- qread(here::here("data/xgboost/finished_model/wrong_vector.qs"))

# saving objects ----
qsave(wang_beta_pre, file = here::here("data/seurat_objects/wang_beta_pre.qs"))

# starting off ---------------------------------------------------------------
# working with things from the script "xgboost_wang"

## making pre seurat object ------------------------------------------
wang_beta_pre <- subset(wang_beta, subset = disease == "pre")
View(wang_beta_pre@meta.data)

## adding prediction to meta data ------------------------------------
# nu vil jeg merge predictions fra M_new (for prediabetics) over i metadata
# alt dette er fra loopet, jeg har genbrugt
temp_label=rep('no',dim(M_new)[1])
temp_pro=rep(-1,dim(M_new)[1])

pred_label=as.numeric(pred > 0.5)
pred_label1=pred_label
pred_label1[pred_label==0]='nd' 
pred_label1[pred_label==1]='t2d' 

keep_test=(as.character(M_new$disease)== "pre")

temp_label[keep_test]=pred_label1
temp_pro[keep_test]=pred

M_new_pre <- M_new # denne er ny, for ikke at ødelægge originalen

M_new_pre$subtype=temp_label 
M_new_pre$subtype_probability=temp_pro
View(M_new_pre)

# nu skal M_new_pre sorteres så det kun er alle pre diabetics
M_new_pre <- subset(M_new_pre, disease == "pre")

# nu har jeg en dataframe som jeg kan bruge til mit meta data i mit pre seurat objekt
wang_beta_pre@meta.data <-  wang_beta_pre@meta.data %>% 
  tibble::rownames_to_column("barcode") %>% 
  dplyr::left_join(y = M_new_pre %>% select(barcode, subtype, subtype_probability),
                   by = "barcode") %>% 
  tibble::column_to_rownames("barcode")


# plot of subtypes distribution -------------------------------------------

table(wang_beta_pre@meta.data$subtype)

ggplot(wang_beta_pre@meta.data, aes(x = orig.ident, fill = subtype)) +
  geom_bar(position = "fill") +
  labs(
    title = "Subtype distribution in prediabetics",
    x = "Donor",
    y = "Percentage"
  ) +
  theme(axis.text.x = element_text(angle = 45))

# clustering of subtypes --------------------------------------------------
DimPlot(wang_beta_pre, reduction = "umap", 
        label = TRUE, group.by = "subtype")

# finding marker genes ----------------------------------------------------
wang_beta_pre <- FindVariableFeatures(wang_beta_pre,
                                  selection.method = "vst", nfeatures = 100)

## laver de 5000 variable gener om til vektorer
var_genes_50 <- VariableFeatures(wang_beta_pre)
var_genes_100 <- VariableFeatures(wang_beta_pre)


DotPlot(
  wang_beta_pre,
  features = var_genes_100, group.by = "subtype"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("Markergenes in subtypes (100 variable features)") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5)
  )


