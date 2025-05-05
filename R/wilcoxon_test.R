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


# Wilcoxon ----------------------------------------------------------------

# wang --------------------------------------------------------------------

#loader data ind (lavet i "wang_diff_genes_and_GOterm.R")
norm_counts_t2d <- qread(here::here("data/goterm/wang/norm_counts_t2d.qs"))

long <- t(norm_counts_t2d %>% dplyr::filter(gene == "INS"))
long <- as.data.frame(long)
colnames(long) = NULL
# lav øverste række til kolonne navne

# assumptions and preliminary tests

TSSEnrichment_long <- isa_dataframe_merged %>%
  dplyr::select(donor_id, method, TSSEnrichment) %>% # vælger kolonner jeg skal arbejde med
  group_by(donor_id, method) %>% # grupperer efter metode og id (for hver donor per metode)
  dplyr::summarise(mean_TSSEnrichment = mean(TSSEnrichment)) %>%  # beregner median for hver donor per metode
  dplyr::ungroup()

TSSEnrichment_wide <- TSSEnrichment_long %>%
  tidyr::pivot_wider(id_cols = donor_id,
                     names_from = method,
                     values_from = mean_TSSEnrichment) # opdeler dataframe i flere koloner så det er 20 rækker og mono+multi for hver

# Summary statistics

TSSEnrichment_long %>%
  group_by(method) %>%
  get_summary_stats(mean_TSSEnrichment, type = "median_iqr")

# Visualization

box_comparison_TSS <- ggpaired(TSSEnrichment_long, x = "method", y = "mean_TSSEnrichment",
                               order = c("mono", "multi"),
                               ylab = "mean_TSSEnrichment", xlab = "method")

box_comparison_TSS

# Assumptions and preliminary tests

TSSEnrichment_wide$differences <- c(TSSEnrichment_wide$mono - TSSEnrichment_wide$multi)

wilcox_assumption_TSSEnrichment <- gghistogram(TSSEnrichment_wide, x = "differences", y = "..density..",
                                               fill = "steelblue",bins = 5, add_density = TRUE)

wilcox_assumption_TSSEnrichment


