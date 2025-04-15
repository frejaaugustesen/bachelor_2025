# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))

set.seed(100)

# saving objects ----
qsave(wang_beta, file = here::here("data/seurat_objects/wang_beta.qs"))
qsave(motakis_beta, file = here::here("data/seurat_objects/motakis_beta.qs"))

qsave(wang_beta_sub, file = here::here("data/seurat_objects/wang_beta_sub.qs"))

qsave(wang_beta_train, file = here::here("data/seurat_objects/wang_beta_train.qs"))
qsave(wang_beta_test, file = here::here("data/seurat_objects/wang_beta_test.qs"))


# load --------------------------------------------------------------------

wang <- qread("/work/bachelor_2025/data/seurat_objects/wang_seurat_QC.qs")
motakis <- qread("/work/bachelor_2025/data/seurat_objects/motakis/islet_all_subset.qs")
wang_beta <- qread("/work/bachelor_2025/data/seurat_objects/wang_beta.qs")

# preprocess --------------------------------------------------------------

head(wang@meta.data)
head(motakis@meta.data)


# to do -------------------------------------------------------------------

# tilføj annoteringer til wang seurat objekt

# subset ----
motakis_beta <- subset(motakis, subset = motakis_anno == "beta")
wang_beta <- subset(wang, subset = wang_anno == "beta")

# Normalize ----
motakis_beta <- NormalizeData(motakis_beta)
wang_beta <- NormalizeData(wang_beta)

# samme gener ----

## finder 5000 variable gener for begge
motakis_beta <- FindVariableFeatures(motakis_beta,
                                selection.method = "vst", nfeatures = 5000)

wang_beta <- FindVariableFeatures(wang_beta,
                                     selection.method = "vst", nfeatures = 5000)

## laver de 5000 variable gener om til vektorer
var_genes_motakis <- VariableFeatures(motakis_beta)
var_genes_wang <- VariableFeatures(wang_beta)

## scale 
motakis_beta <- ScaleData(motakis_beta, features = var_genes_motakis)
wang_beta <- ScaleData(wang_beta, features = var_genes_wang)


## sammenligner hvilke gener som er i begge vektorer så man beholder kun dem der er begge steder
var_genes <- intersect(var_genes_wang, var_genes_motakis) # kun beholdt 886 gener?

# add disease to wang meta data ----------------------------------------------

meta_wang <- read.table(here::here("data_raw/wang/meta_donor.csv"), 
                   header = TRUE, 
                   sep = "", 
                   stringsAsFactors = FALSE)

# seeing if i need all of the donors in the meta data or if some needs removing
length(intersect(wang_beta@meta.data$orig.ident, meta$dono))
  # all donors in meta data is in the seurat object

# adding disease state to meta data in seurat

#renaming column to be the same in both dataframes (to use left_jointly)
meta_wang <- dplyr::rename(meta_wang, orig.ident = donor)

# add donor meta data to seurat object
wang_beta@meta.data <-  wang_beta@meta.data %>% 
  tibble::rownames_to_column("barcode") %>% 
  dplyr::left_join(y = meta_wang, by = "orig.ident") %>% 
  tibble::column_to_rownames("barcode")

# first attempt ---------------------------------------------------------------

## extract pre diabetics ------------------------------------------------------
wang_beta_sub <- subset(wang_beta, subset = disease != "pre")

# checking they were removed (before and after)
table(wang_beta@meta.data$disease)
table(wang_beta_sub@meta.data$disease)

## divide data into model data and data for subtype prediction ----------------
length(table(wang_beta_sub@meta.data$orig.ident)) # der er 12 donor nu

# getting all cells
all_cells <- colnames(wang_beta_sub)

# randomly sample 80% 
train_cells <- sample(all_cells, size = 0.8 * length(all_cells))

# remaining 20% go to test set
test_cells <- setdiff(all_cells, train_cells)

# subset the Seurat objects intro train and test
wang_beta_train <- subset(wang_beta_sub, cells = train_cells)
wang_beta_test <- subset(wang_beta_sub, cells = test_cells)

## make structure for model ---------------------------------------------------

### data ----
# need to be a dgCMatrix or a matrix, with variables and the cells

#### train ----
# using only var_genes
# Get the log-normalized data matrix 
rna_data_train <- GetAssayData(wang_beta_train, assay = "RNA", slot = "data")

# keeping genes, which are in var_genes (determined earlier)
matrix_data_train <- rna_data_train[var_genes, ]

length(matrix_data_train)

#### test ----
# using only var_genes
# Get the log-normalized data matrix 
rna_data_test <- GetAssayData(wang_beta_test, assay = "RNA", slot = "data")

# keeping genes, which are in var_genes (determined earlier)
matrix_data_test <- rna_data_test[var_genes, ]

length(matrix_data_test)

# now the matrix for the data is created for train and test
matrix_data_train
matrix_data_test

### label ----
# needs to bee of class "numeric", as 0 and 1 for disease state
disease_vector_test <- wang_beta_test@meta.data$disease
disease_vector_train <- wang_beta_train@meta.data$disease

label_test <- recode(disease_vector_test, "nd" = 0, "t2d" = 1)
label_train <- recode(disease_vector_train, "nd" = 0, "t2d" = 1)

### combine ----
# now i have the data and the label, which need to combined in a list
train_wang <- list(matrix = matrix_data_train, vector = label_train)
test_wang <- list(matrix = matrix_data_test, vector = label_test)

## conclusion ----------------------------------------------------------------
# now the structure of the data is like in the mushroom vignette, but donor id's 
# can not be seen in the data, and therefore not run in a loop and perform
# leave one out analyzis. 


# second attempt ----------------------------------------------------------
# i will try to attempt to recreate code from the wang sander study
# https://github.com/gaoweiwang/Islet_snATACseq/blob/8145d545d979507193fd5bef71b1b4119bed6e66/scripts/AI_subtype.ipynb#L710

## setting up data ----
# first creating dataframe with disease and predicted disease 
M_new <- wang_beta@meta.data
M_new$pre_disease <- wang_beta@meta.data$disease
M_new <- dplyr::rename(M_new, donor = orig.ident)


# creating vector with all donor names
wang_beta_sub <- subset(wang_beta, subset = disease != "pre") # used earlier
wang_meta_T2D_ND <- wang_beta_sub@meta.data
donor_all <- unique(wang_meta_T2D_ND$orig.ident)

# hvad betyder det her??
temp_label=rep('no',dim(M_new)[1])
temp_pro=rep(-1,dim(M_new)[1])

## setting up loop ----
# arbejder på at forstå dette
for (i in 1:length(donor_all)){
  keep_1=(as.character(M_new$donor)!=donor_all[i])
  keep_2=(as.character(M_new$disease)!='pd')
  keep_3=(as.character(M_new$disease)==as.character(M_new$pre_disease))
  
  keep_test=(as.character(M_new$donor)==donor_all[i])
  
  keep_train=(keep_1 & keep_2 & keep_3)
  train.x = data_use[keep_train,]
  test.x = data_use[keep_test,]
  
  train_disease=as.character(M_new$disease)
  train_d=rep(0,length(train_disease))
  train_d[train_disease=='t2d']=1   #######
  train_d[train_disease=='nd']=0  #######
  train.y = train_d[keep_train]
  test.y = train_d[keep_test]
  
  train.x1=as(train.x, "dgCMatrix")
  bst <- xgboost(data = train.x1, label = train.y, max.depth = 60, eta = 0.2, nthread = 24, nrounds = 80, objective = "binary:logistic")
  test.x1=as(test.x, "dgCMatrix")
  pred <- predict(bst, test.x1)
  
  pred_label=as.numeric(pred > 0.5)
  pred_label1=pred_label
  pred_label1[pred_label==0]='nd' #######
  pred_label1[pred_label==1]='t2d' #######
  temp_label[keep_test]=pred_label1
  temp_pro[keep_test]=pred
  
  M_new$subtype=temp_label #######
  M_new$subtype_probability=temp_pro
  write.csv(M_new,paste0(wd_snATAC,'beta_subtype_temp2.csv')) #######
}





