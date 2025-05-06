# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group
# testing for significance in subtypes and gene expression

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))

set.seed(100)


# load --------------------------------------------------------------------

motakis_beta <- qread("/work/bachelor_2025/data/seurat_objects/motakis_beta_integrated.qs")
wang_beta <- qread("/work/bachelor_2025/data/seurat_objects/wang_beta_pred.qs")


# Chi-squared test --------------------------------------------------------

## Wang --------------------------------------------------------------------

table(wang_beta@meta.data$subtype, wang_beta@meta.data$disease)
# hvad skal jeg sammenligne helt præcist??

# Create the table
wang_table <- table(wang_beta@meta.data$subtype, wang_beta@meta.data$disease)

# Chi-squared test
chisq.test(wang_table)
# der er en sammenhæng mellem subtype og disease fordi p-værdi er lav


## Motakis -----------------------------------------------------------------

table(motakis_beta@meta.data$subtype, motakis_beta@meta.data$disease)
# hvad skal jeg sammenligne helt præcist??

# Create the table
motakis_table <- table(motakis_beta@meta.data$subtype, motakis_beta@meta.data$disease)

# Chi-squared test
chisq.test(motakis_table)
# der er en sammenhæng mellem subtype og disease fordi p-værdi er lav


