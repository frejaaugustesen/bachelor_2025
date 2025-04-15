# omdøb filen Freja!
# arbejder med hvordan man samler seuratobjekter og harmony
merged_seurat <- merge(islet56, y= islet57, add.cell.ids = c("islet56", "islet57"))

str(merged_seurat)

# harmony
library(harmony)

# man merger alle og så splitter sit RNA lag ud fra en faktor (f.eks. orig.ident - så den splitter per prøve)
islet56 <- Seurat::CreateSeuratObject(counts = mtx_gene_sub,
                                      assay = "RNA", project = "islet56")
# her bruges project til at kalde orig.ident for islet56 istedet for seuratproject

# i eksempel bruger de integrated.caa - her skal jeg bruge harmony.