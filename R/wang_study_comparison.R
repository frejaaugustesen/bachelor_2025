# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group
# comparing with wang study

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))

set.seed(100)


# load --------------------------------------------------------------------

wang_beta <- qread("/work/bachelor_2025/data/seurat_objects/wang_beta_pred.qs")
study <- read.csv(here::here("data_raw/wang/multiome_RNA_beta_subtype.csv"))


# renaming subtypes -------------------------------------------------------

#renaming subtypes in meta data from nd and t2d to beta1 and beta2
wang_beta@meta.data <- wang_beta@meta.data %>%
  mutate(subtype = recode(subtype,
                          "nd" = "beta1",
                          "t2d" = "beta2"))

# adding study subtypes to metadata ---------------------------------------------

# column to rownames in study csv

study$X <- sub("-1$", "", study$X)

study <- study %>%
  rename(barcode = X) %>%
  rename(study_subtype = subtype)


# adding the wang study subtype to meta data
wang_beta@meta.data <-  wang_beta@meta.data %>% 
  tibble::rownames_to_column("barcode") %>% 
  dplyr::left_join(y = study %>% select(barcode, study_subtype),
                   by = "barcode") %>% 
  mutate(agreement = case_when(subtype == study_subtype ~ "yes", 
                               subtype != study_subtype ~ "no")) %>% 
  tibble::column_to_rownames("barcode")


# comparing ---------------------------------------------------------------

table(wang_beta@meta.data$agreement)

DimPlot(wang_beta, reduction = "umap", label = TRUE,
        group.by = c("subtype", "study_subtype", "agreement"),
        label.size = 2.5,
        repel = TRUE) &
  NoLegend()

DimPlot(wang_beta, reduction = "umap", label = TRUE,
        group.by = "agreement",
        label.size = 2.5,
        repel = TRUE)



