# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))

set.seed(100)


# load --------------------------------------------------------------------

motakis_beta <- qread("/work/bachelor_2025/data/seurat_objects/motakis_beta_integrated.qs")
wang_beta <- qread("/work/bachelor_2025/data/seurat_objects/wang_beta_pred.qs")

meta_data <- read.csv("/work/bachelor_2025/data_raw/motakis/motakis_meta.csv", sep = ";") %>% 
  select(-disease)

wang_meta <- read.csv("/work/bachelor_2025/data_raw/wang/meta_donor.csv", sep = "\t")

# laver om til procent
motakis_tabel <- motakis_beta@meta.data %>% 
  as.data.frame() %>% 
  group_by(disease, subtype, orig.ident, sample) %>% 
  tally() %>% 
  group_by(orig.ident) %>% 
  mutate(perc = (n/sum(n)*100), 
         sum = sum(n)) %>% 
  filter(subtype == "t2d") %>% 
  left_join(y = meta_data, by = "sample")
# bruges til at lave korrelations plot

ggscatter(motakis_tabel, x = "perc", y = "hba1c_.",
          add = "reg.line") +
  stat_cor(method = "pearson")

ggscatter(motakis_tabel, x = "perc", y = "bmi",
          add = "reg.line") +
  stat_cor(method = "pearson")

ggscatter(motakis_tabel, x = "perc", y = "age_years",
          add = "reg.line") +
  stat_cor(method = "pearson")


## gentag for wang
# laver om til procent
wang_tabel <- wang_beta@meta.data %>% 
  as.data.frame() %>% 
  rename(donor = orig.ident) %>% 
  group_by(disease, subtype, donor) %>% 
  tally() %>% 
  group_by(donor) %>% 
  mutate(perc = (n/sum(n)*100), 
         sum = sum(n)) %>% 
  filter(subtype == "t2d") %>% 
  left_join(y = wang_meta, by = "donor")
# bruges til at lave korrelations plot

ggscatter(wang_tabel, x = "perc", y = "hba1c_.",
          add = "reg.line") +
  stat_cor(method = "pearson")

ggscatter(wang_tabel, x = "perc", y = "bmi",
          add = "reg.line") +
  stat_cor(method = "pearson")

ggscatter(wang_tabel, x = "perc", y = "age_years",
          add = "reg.line") +
  stat_cor(method = "pearson") 
# find ud af at farve efter disease
  scale_fill_manual(values = c(nd = "pink", pre = "lightgreen", t2d = "lightblue"))
  

# gener -------------------------------------------------------------------

motakis_norm_counts <- norm_counts #defineret i "motakis_diff_genes_and_GOterm"
wang_res_sig <- res_sig ##defineret i "wang_diff_genes_and_GOterm"

wang_res_sig <- res %>%
  as.data.frame() %>%
  dplyr::filter(padj <= 0.05 & abs(log2FoldChange) > 1) %>% 
  rownames()


motakis_genes <- filter(motakis_norm_counts, gene %in% wang_res_sig) %>% 
  tibble::column_to_rownames("gene") # brug denne på wang

motakis_genes %>% 
  rowwise() %>% 
  filter(colnames(motakis_genes) != 0)

motakis_heat <- motakis_genes %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum != 0) %>% 
  select(-sum)

motakis_genes$sum

pheatmap::pheatmap(motakis_heat, 
                   scale = "row",
                   color = myCol,
                   breaks = myBreaks,
                   cluster_cols = FALSE) # brug denne på wang


norm_counts %>% dplyr::filter(gene == "NTNG1") %>%
  tidyr::pivot_longer(-gene, names_to = "sample", values_to = "exp") %>%
  dplyr::mutate(subtype = case_when(grepl("nd", sample) ~ "nd",
                                    grepl("t2d", sample) ~ "t2d")) %>%
  ggplot2::ggplot(aes(x = subtype, y = exp)) +
  geom_bar(stat='summary', fun = "mean") +
  ggplot2::geom_point()+
  ggtitle("Motakis subtypes, gene NTNG1")

# wang

wang_norm_counts <- norm_counts
wang_genes <- filter(wang_norm_counts, gene %in% wang_res_sig) %>% 
  tibble::column_to_rownames("gene") # brug denne på wang

myCol <- colorRampPalette(c('#004B7A', 'white', '#A83708'))(100)
myBreaks <- seq(-1.5, 1.5, length.out = 100)

pheatmap::pheatmap(wang_genes, 
                   scale = "row",
                   color = myCol,
                   breaks = myBreaks) # brug denne på wang
