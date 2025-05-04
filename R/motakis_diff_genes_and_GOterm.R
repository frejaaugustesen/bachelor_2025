# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group
# differentially expressed genes and GO-term analzysis for motakis

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))

set.seed(100)

# load --------------------------------------------------------------------

motakis_beta <- qread("/work/bachelor_2025/data/seurat_objects/motakis_beta_integrated.qs")

# Differentially expressed genes ------------------------------------------
# Install needed packages
BiocManager::install("edgeR")
BiocManager::install("DESeq2")

# Subset to only contain real T2D individuals 
# (Jeg tror måske det er bedre at sammenligne subtyper, for hver disease?)
Idents(motakis_beta) <- "disease"
motakis_beta_sub <- subset(x = motakis_beta, idents = "t2d")

# Prepare data for differential gene expression analysis
counts <- motakis_beta_sub %>%
  edgeR::Seurat2PB(sample = "orig.ident",
                   cluster = "subtype")

# Extract pseudobulk counts
counts_2 <- counts[["counts"]] %>%
  as.data.frame() %>%
  dplyr::rename_with(~gsub("cluster", "", .x))

# Get meta data
meta_data <- counts[["samples"]] %>%
  magrittr::set_rownames(base::gsub("cluster", "", BiocGenerics::rownames(.)))

# Check colnames and rownames are equal
base::all.equal(BiocGenerics::colnames(counts_2), BiocGenerics::rownames(meta_data))

# Create a deseq2 object - paired test (this is why we include sample)
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_2,
  colData = meta_data,
  design = stats::as.formula("~ sample + cluster")
)

# Normalize and run DESeq
dds <- DESeq2::DESeq(dds)

# Find diff genes between nd and t2d subtypes (gene uprgulatedin nd vs t2d)
res <- DESeq2::results(dds, contrast = c("cluster", "nd", "t2d"))

# Significant results
res_sig <- res %>%
  as.data.frame() %>%
  dplyr::filter(padj <= 0.05)

# We do not see any significant genes but that also makes sense if we look at the PCA plot:

# Transform using regularised logarithm
dds_rlog <- DESeq2::rlogTransformation(dds)

# Create PCA plot
DESeq2::plotPCA(dds_rlog, intgroup="cluster")


# Get normalized counts ---------------------------------------------------
norm_counts <- DESeq2::counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene")

# Plot expression of your favorite gene
norm_counts %>% dplyr::filter(gene == "INS") %>%
  tidyr::pivot_longer(-gene, names_to = "sample", values_to = "exp") %>%
  dplyr::mutate(subtype = case_when(grepl("nd", sample) ~ "nd",
                                    grepl("t2d", sample) ~ "t2d")) %>%
  ggplot2::ggplot(aes(x = subtype, y = exp)) +
  geom_bar(stat='summary', fun = "mean") +
  ggplot2::geom_point()

# Go_term_analysis --------------------------------------------------------
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

# we just remove genes where the p-value is NA, although we should not do this
# for the real analysis (but we dont have any significant ones)
res_sig <- res %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(padj))

# Look both at up and downregulated genes separately - this example is for
# upregulated genes
# (it could also be above 0) - these are genes upregulated in nd subtypes
res_up <- res_sig %>%
  dplyr::filter(log2FoldChange > 1)

# Convert gene symbols to entrez ids
res_up$entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                      keys=rownames(res_up),
                                      column="ENTREZID",
                                      keytype="SYMBOL",
                                      multiVals="first")

# Background genes - all genes that were tested
res$entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                   keys=rownames(res),
                                   column="ENTREZID",
                                   keytype="SYMBOL",
                                   multiVals="first")
# Universe
universe <- res %>%
  as.data.frame() %>%
  filter(!is.na(entrez)) %>%
  pull(entrez) %>%
  unlist() %>%
  unname() %>%
  unique()

# genes
genes <- res_up %>%
  as.data.frame() %>%
  filter(!is.na(entrez)) %>%
  pull(entrez) %>%
  unlist() %>%
  unname() %>%
  unique()

# go-term for biological ontologies
go_up <- clusterProfiler::enrichGO(
  gene = genes,
  OrgDb = "org.Hs.eg.db",
  pAdjustMethod = "fdr",
  ont = "BP",
  pvalueCutoff  = 0.2,
  qvalueCutoff  = 0.2,
  readable      = TRUE,
  universe = universe)

# Get top 5 go-term
go_up@result %>%
  dplyr::filter(p.adjust <= 0.05) %>%
  dplyr::arrange(p.adjust) %>%
  head(n = 5)
