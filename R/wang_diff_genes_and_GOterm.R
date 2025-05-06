# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group
# differentially expressed genes and GO-term analzysis for wang

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))

set.seed(100)

# load --------------------------------------------------------------------

wang_beta <- qread("/work/bachelor_2025/data/seurat_objects/wang_beta_pred.qs")


# save --------------------------------------------------------------------

#qsave(top5_up_t2d, file = here::here("data/goterm/top5_up_t2d.qs"))
#qsave(top5_down_t2d, file = here::here("data/goterm/top5_down_t2d.qs"))

#qsave(top5_up_nd, file = here::here("data/goterm/top5_up_nd.qs"))
#qsave(top5_down_nd, file = here::here("data/goterm/top5_down_nd.qs"))

#qsave(top5_up_pre, file = here::here("data/goterm/top5_up_pre.qs"))
#qsave(top5_down_pre, file = here::here("data/goterm/top5_down_pre.qs"))

# T2D donors --------------------------------------------------------------

# Subset to only contain real T2D individuals 
# (Jeg tror måske det er bedre at sammenligne subtyper, for hver disease?)
Idents(wang_beta) <- "disease"
wang_beta_t2d <- subset(x = wang_beta, idents = "t2d")


## Marker genes ------------------------------------------------------------
Idents(wang_beta_t2d) <- "subtype"

markergenes_t2d <- FindMarkers(wang_beta_t2d, ident.1 = "nd", ident.2 = "t2d", 
                           group.by = "subtype")



allmarkers_t2d <- FindAllMarkers(wang_beta_t2d, only.pos = TRUE, min.pct = 0.1)


markergenes_t2d_filt <- markergenes_t2d %>%
  filter(p_val_adj <= 0.05, avg_log2FC > 0) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  tibble::rownames_to_column("gene")

allmarkers_t2d_filt <- allmarkers_t2d %>%
  filter(p_val_adj <= 0.05) %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)


DotPlot(
  wang_beta_t2d,
  features = allmarkers_t2d_filt$gene, group.by = "subtype"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("Markergenes in T2D subtypes (FindAllMarkers,only.pos = TRUE, min.pct = 0.1)") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5)
  )

## Differentially expressed genes ------------------------------------------

Idents(wang_beta_t2d) <- "disease"

# Prepare data for differential gene expression analysis
counts_t2d <- wang_beta_t2d %>%
  edgeR::Seurat2PB(sample = "orig.ident",
                   cluster = "subtype")

# Extract pseudobulk counts
counts_t2d_2 <- counts_t2d[["counts"]] %>%
  as.data.frame() %>%
  dplyr::rename_with(~gsub("cluster", "", .x))

# Get meta data
meta_data_t2d <- counts_t2d[["samples"]] %>%
  magrittr::set_rownames(base::gsub("cluster", "", BiocGenerics::rownames(.)))

# Check colnames and rownames are equal
base::all.equal(BiocGenerics::colnames(counts_t2d_2), BiocGenerics::rownames(meta_data_t2d))

# Create a deseq2 object - paired test (this is why we include sample)
dds_t2d <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_t2d_2,
  colData = meta_data_t2d,
  design = stats::as.formula("~ sample + cluster")
)

# Normalize and run DESeq
dds_t2d <- DESeq2::DESeq(dds_t2d)

# Find diff genes between nd and t2d subtypes (gene uprgulatedin nd vs t2d)
res_t2d <- DESeq2::results(dds_t2d, contrast = c("cluster", "nd", "t2d"))

# Significant results
res_sig_t2d <- res_t2d %>%
  as.data.frame() %>%
  dplyr::filter(padj <= 0.05)

# Transform using regularised logarithm
dds_rlog_t2d <- DESeq2::rlogTransformation(dds_t2d)

# Create PCA plot
DESeq2::plotPCA(dds_rlog_t2d, intgroup="cluster")


## Get normalized counts ---------------------------------------------------
norm_counts_t2d <- DESeq2::counts(dds_t2d, normalized = TRUE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene")

# Plot expression of your favorite gene
norm_counts_t2d %>% dplyr::filter(gene == "ROR1") %>%
  tidyr::pivot_longer(-gene, names_to = "sample", values_to = "exp") %>%
  dplyr::mutate(subtype = case_when(grepl("nd", sample) ~ "nd",
                                    grepl("t2d", sample) ~ "t2d")) %>%
  ggplot2::ggplot(aes(x = subtype, y = exp)) +
  geom_bar(stat='summary', fun = "mean") +
  ggplot2::geom_point() +
  labs(title = "Expression of INS in subtypes")

qsave(norm_counts_t2d, here::here("data/goterm/wang/norm_counts_t2d.qs"))

## Go_term_analysis --------------------------------------------------------

# Look both at up and downregulated genes separately - this example is for
# upregulated genes
# (it could also be above 0) - these are genes upregulated in nd subtypes
res_up_t2d <- res_sig_t2d %>%
  dplyr::filter(log2FoldChange > 0)

res_down_t2d <- res_sig_t2d %>%
  dplyr::filter(log2FoldChange < 0)

# Convert gene symbols to entrez ids
res_up_t2d$entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                      keys=rownames(res_up_t2d),
                                      column="ENTREZID",
                                      keytype="SYMBOL",
                                      multiVals="first")

res_down_t2d$entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                          keys=rownames(res_down_t2d),
                                          column="ENTREZID",
                                          keytype="SYMBOL",
                                          multiVals="first")

# Background genes - all genes that were tested
res_t2d$entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                   keys=rownames(res_t2d),
                                   column="ENTREZID",
                                   keytype="SYMBOL",
                                   multiVals="first")
# Universe
universe_t2d <- res_t2d %>%
  as.data.frame() %>%
  filter(!is.na(entrez)) %>%
  pull(entrez) %>%
  unlist() %>%
  unname() %>%
  unique()

# genes up
genes_up_t2d <- res_up_t2d %>%
  as.data.frame() %>%
  filter(!is.na(entrez)) %>%
  pull(entrez) %>%
  unlist() %>%
  unname() %>%
  unique()

# genes down
genes_down_t2d <- res_down_t2d %>%
  as.data.frame() %>%
  filter(!is.na(entrez)) %>%
  pull(entrez) %>%
  unlist() %>%
  unname() %>%
  unique()

# go-term for biological ontologies
go_up_t2d <- clusterProfiler::enrichGO(
  gene = genes_up_t2d,
  OrgDb = "org.Hs.eg.db",
  pAdjustMethod = "fdr",
  ont = "BP",
  pvalueCutoff  = 0.2,
  qvalueCutoff  = 0.2,
  readable      = TRUE,
  universe = universe_t2d)

# go-term for biological ontologies
go_down_t2d <- clusterProfiler::enrichGO(
  gene = genes_down_t2d,
  OrgDb = "org.Hs.eg.db",
  pAdjustMethod = "fdr",
  ont = "BP",
  pvalueCutoff  = 0.2,
  qvalueCutoff  = 0.2,
  readable      = TRUE,
  universe = universe_t2d)

# Get top 5 go-term
top5 <- go_up_t2d@result %>%
  dplyr::filter(p.adjust <= 0.05) %>% # hvorfor skal den være under 0.05??
  dplyr::arrange(p.adjust) %>%
  head(n = 5)

top5 <- go_down_t2d@result %>%
  dplyr::filter(p.adjust <= 0.05) %>% # hvorfor skal den være under 0.05??
  dplyr::arrange(p.adjust) %>%
  head(n = 5)

top5_up_t2d <- go_up_t2d@result %>%
  dplyr::arrange(p.adjust) %>%
  head(n = 5)

top5_down_t2d <- go_down_t2d@result %>%
  dplyr::arrange(p.adjust) %>%
  head(n = 5)

# ND donors --------------------------------------------------------------

# Subset to only contain real T2D individuals 
# (Jeg tror måske det er bedre at sammenligne subtyper, for hver disease?)
wang_beta_nd <- subset(x = wang_beta, idents = "nd")


## Marker genes ------------------------------------------------------------
Idents(wang_beta_nd) <- "subtype"

markergenes_nd <- FindMarkers(wang_beta_nd, ident.1 = "nd", ident.2 = "t2d", 
                               group.by = "subtype")



allmarkers_nd <- FindAllMarkers(wang_beta_nd, only.pos = TRUE, min.pct = 0.1)


markergenes_nd_filt <- markergenes_nd %>%
  filter(p_val_adj <= 0.05, avg_log2FC > 0) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  tibble::rownames_to_column("gene")

allmarkers_nd_filt <- allmarkers_nd %>%
  filter(p_val_adj <= 0.05) %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)


DotPlot(
  wang_beta_nd,
  features = allmarkers_nd_filt$gene, group.by = "subtype"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("Markergenes in ND subtypes (FindAllMarkers,only.pos = TRUE, min.pct = 0.1)") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5)
  )

## Differentially expressed genes ------------------------------------------

Idents(wang_beta_nd) <- "disease"

# Prepare data for differential gene expression analysis
counts_nd <- wang_beta_nd %>%
  edgeR::Seurat2PB(sample = "orig.ident",
                   cluster = "subtype")

# Extract pseudobulk counts
counts_nd_2 <- counts_nd[["counts"]] %>%
  as.data.frame() %>%
  dplyr::rename_with(~gsub("cluster", "", .x))

# Get meta data
meta_data_nd <- counts_nd[["samples"]] %>%
  magrittr::set_rownames(base::gsub("cluster", "", BiocGenerics::rownames(.)))

# Check colnames and rownames are equal
base::all.equal(BiocGenerics::colnames(counts_nd_2), BiocGenerics::rownames(meta_data_nd))

# Create a deseq2 object - paired test (this is why we include sample)
dds_nd <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_nd_2,
  colData = meta_data_nd,
  design = stats::as.formula("~ sample + cluster")
)

# Normalize and run DESeq
dds_nd <- DESeq2::DESeq(dds_nd)

# Find diff genes between nd and t2d subtypes (gene uprgulatedin nd vs t2d)
res_nd <- DESeq2::results(dds_nd, contrast = c("cluster", "nd", "t2d"))

# Significant results
res_sig_nd <- res_nd %>%
  as.data.frame() %>%
  dplyr::filter(padj <= 0.05)

# Transform using regularised logarithm
dds_rlog_nd <- DESeq2::rlogTransformation(dds_nd)

# Create PCA plot
DESeq2::plotPCA(dds_rlog_nd, intgroup="cluster")


## Get normalized counts ---------------------------------------------------
norm_counts_nd <- DESeq2::counts(dds_nd, normalized = TRUE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene")

# Plot expression of your favorite gene
norm_counts_nd %>% dplyr::filter(gene == "INS") %>%
  tidyr::pivot_longer(-gene, names_to = "sample", values_to = "exp") %>%
  dplyr::mutate(subtype = case_when(grepl("nd", sample) ~ "nd",
                                    grepl("t2d", sample) ~ "t2d")) %>%
  ggplot2::ggplot(aes(x = subtype, y = exp)) +
  geom_bar(stat='summary', fun = "mean") +
  ggplot2::geom_point() +
  labs(title = "Expression of INS in subtypes")

## Go_term_analysis --------------------------------------------------------

# Look both at up and downregulated genes separately - this example is for
# upregulated genes
# (it could also be above 0) - these are genes upregulated in nd subtypes
res_up_nd <- res_sig_nd %>%
  dplyr::filter(log2FoldChange > 0.0)

res_down_nd <- res_sig_nd %>%
  dplyr::filter(log2FoldChange < 0.0)

# Convert gene symbols to entrez ids
res_up_nd$entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                          keys=rownames(res_up_nd),
                                          column="ENTREZID",
                                          keytype="SYMBOL",
                                          multiVals="first")

res_down_nd$entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                            keys=rownames(res_down_nd),
                                            column="ENTREZID",
                                            keytype="SYMBOL",
                                            multiVals="first")

# Background genes - all genes that were tested
res_nd$entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                       keys=rownames(res_nd),
                                       column="ENTREZID",
                                       keytype="SYMBOL",
                                       multiVals="first")
# Universe
universe_nd <- res_nd %>%
  as.data.frame() %>%
  filter(!is.na(entrez)) %>%
  pull(entrez) %>%
  unlist() %>%
  unname() %>%
  unique()

# genes up
genes_up_nd <- res_up_nd %>%
  as.data.frame() %>%
  filter(!is.na(entrez)) %>%
  pull(entrez) %>%
  unlist() %>%
  unname() %>%
  unique()

# genes down
genes_down_nd <- res_down_nd %>%
  as.data.frame() %>%
  filter(!is.na(entrez)) %>%
  pull(entrez) %>%
  unlist() %>%
  unname() %>%
  unique()

# go-term for biological ontologies
go_up_nd <- clusterProfiler::enrichGO(
  gene = genes_up_nd,
  OrgDb = "org.Hs.eg.db",
  pAdjustMethod = "fdr",
  ont = "BP",
  pvalueCutoff  = 0.2,
  qvalueCutoff  = 0.2,
  readable      = TRUE,
  universe = universe_nd)

# go-term for biological ontologies
go_down_nd <- clusterProfiler::enrichGO(
  gene = genes_down_nd,
  OrgDb = "org.Hs.eg.db",
  pAdjustMethod = "fdr",
  ont = "BP",
  pvalueCutoff  = 0.2,
  qvalueCutoff  = 0.2,
  readable      = TRUE,
  universe = universe_nd)

# Get top 5 go-term
top5 <- go_up_nd@result %>%
  dplyr::filter(p.adjust <= 0.05) %>% # hvorfor skal den være under 0.05??
  dplyr::arrange(p.adjust) %>%
  head(n = 5)

top5_up_nd <- go_up_nd@result %>%
  dplyr::arrange(p.adjust) %>%
  head(n = 5)

top5_down_nd <- go_down_nd@result %>%
  dplyr::arrange(p.adjust) %>%
  head(n = 5)


# PRE donors --------------------------------------------------------------

# Subset to only contain real T2D individuals 
# (Jeg tror måske det er bedre at sammenligne subtyper, for hver disease?)
wang_beta_pre <- subset(x = wang_beta, idents = "pre")


## Marker genes ------------------------------------------------------------
Idents(wang_beta_pre) <- "subtype"

markergenes_pre <- FindMarkers(wang_beta_pre, ident.1 = "nd", ident.2 = "t2d", 
                              group.by = "subtype")



allmarkers_pre <- FindAllMarkers(wang_beta_pre, only.pos = TRUE, min.pct = 0.1)


markergenes_pre_filt <- markergenes_pre %>%
  filter(p_val_adj <= 0.05, avg_log2FC > 0) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  tibble::rownames_to_column("gene")

allmarkers_pre_filt <- allmarkers_pre %>%
  filter(p_val_adj <= 0.05) %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)


DotPlot(
  wang_beta_pre,
  features = allmarkers_pre_filt$gene, group.by = "subtype"
) +
  ggplot2::scale_colour_gradient2(
    low = "#004B7AFF",
    mid = "#FDFDFCFF",
    high = "#A83708FF"
  ) +
  ggtitle("Markergenes in pre subtypes (FindAllMarkers,only.pos = TRUE, min.pct = 0.1)") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5)
  )

## Differentially expressed genes ------------------------------------------

Idents(wang_beta_pre) <- "disease"

# Prepare data for differential gene expression analysis
counts_pre <- wang_beta_pre %>%
  edgeR::Seurat2PB(sample = "orig.ident",
                   cluster = "subtype")

# Extract pseudobulk counts
counts_pre_2 <- counts_pre[["counts"]] %>%
  as.data.frame() %>%
  dplyr::rename_with(~gsub("cluster", "", .x))

# Get meta data
meta_data_pre <- counts_pre[["samples"]] %>%
  magrittr::set_rownames(base::gsub("cluster", "", BiocGenerics::rownames(.)))

# Check colnames and rownames are equal
base::all.equal(BiocGenerics::colnames(counts_pre_2), BiocGenerics::rownames(meta_data_pre))

# Create a deseq2 object - paired test (this is why we include sample)
dds_pre <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_pre_2,
  colData = meta_data_pre,
  design = stats::as.formula("~ sample + cluster")
)

# Normalize and run DESeq
dds_pre <- DESeq2::DESeq(dds_pre)

# Find diff genes between nd and t2d subtypes (gene uprgulatedin nd vs t2d)
res_pre <- DESeq2::results(dds_pre, contrast = c("cluster", "nd", "t2d"))

# Significant results
res_sig_pre <- res_pre %>%
  as.data.frame() %>%
  dplyr::filter(padj <= 0.05)

# Transform using regularised logarithm
dds_rlog_pre <- DESeq2::rlogTransformation(dds_pre)

# Create PCA plot
DESeq2::plotPCA(dds_rlog_pre, intgroup="cluster")


## Get normalized counts ---------------------------------------------------
norm_counts_pre <- DESeq2::counts(dds_pre, normalized = TRUE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene")

# Plot expression of your favorite gene
norm_counts_pre %>% dplyr::filter(gene == "INS") %>%
  tidyr::pivot_longer(-gene, names_to = "sample", values_to = "exp") %>%
  dplyr::mutate(subtype = case_when(grepl("nd", sample) ~ "nd",
                                    grepl("t2d", sample) ~ "t2d")) %>%
  ggplot2::ggplot(aes(x = subtype, y = exp)) +
  geom_bar(stat='summary', fun = "mean") +
  ggplot2::geom_point() +
  labs(title = "Expression of INS in subtypes")

## Go_term_analysis --------------------------------------------------------

# Look both at up and downregulated genes separately - this example is for
# upregulated genes
# (it could also be above 0) - these are genes upregulated in nd subtypes
res_up_pre <- res_sig_pre %>%
  dplyr::filter(log2FoldChange > 1)

res_down_pre <- res_sig_pre %>%
  dplyr::filter(log2FoldChange < 1)

# Convert gene symbols to entrez ids
res_up_pre$entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                         keys=rownames(res_up_pre),
                                         column="ENTREZID",
                                         keytype="SYMBOL",
                                         multiVals="first")

res_down_pre$entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                           keys=rownames(res_down_pre),
                                           column="ENTREZID",
                                           keytype="SYMBOL",
                                           multiVals="first")

# Background genes - all genes that were tested
res_pre$entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                      keys=rownames(res_pre),
                                      column="ENTREZID",
                                      keytype="SYMBOL",
                                      multiVals="first")
# Universe
universe_pre <- res_pre %>%
  as.data.frame() %>%
  filter(!is.na(entrez)) %>%
  pull(entrez) %>%
  unlist() %>%
  unname() %>%
  unique()

# genes up
genes_up_pre <- res_up_pre %>%
  as.data.frame() %>%
  filter(!is.na(entrez)) %>%
  pull(entrez) %>%
  unlist() %>%
  unname() %>%
  unique()

# genes down
genes_down_pre <- res_down_pre %>%
  as.data.frame() %>%
  filter(!is.na(entrez)) %>%
  pull(entrez) %>%
  unlist() %>%
  unname() %>%
  unique()

# go-term for biological ontologies
go_up_pre <- clusterProfiler::enrichGO(
  gene = genes_up_pre,
  OrgDb = "org.Hs.eg.db",
  pAdjustMethod = "fdr",
  ont = "BP",
  pvalueCutoff  = 0.2,
  qvalueCutoff  = 0.2,
  readable      = TRUE,
  universe = universe_pre)

# go-term for biological ontologies
go_down_pre <- clusterProfiler::enrichGO(
  gene = genes_down_pre,
  OrgDb = "org.Hs.eg.db",
  pAdjustMethod = "fdr",
  ont = "BP",
  pvalueCutoff  = 0.2,
  qvalueCutoff  = 0.2,
  readable      = TRUE,
  universe = universe_pre)

# Get top 5 go-term
top5 <- go_up_pre@result %>%
  dplyr::filter(p.adjust <= 0.05) %>% # hvorfor skal den være under 0.05??
  dplyr::arrange(p.adjust) %>%
  head(n = 5)

top5_up_pre <- go_up_pre@result %>%
  dplyr::arrange(p.adjust) %>%
  head(n = 5)

top5_down_pre <- go_down_pre@result %>%
  dplyr::arrange(p.adjust) %>%
  head(n = 5)


