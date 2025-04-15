# genes expressed by different cell types

# pancreas -------------------------------------------------------------------
beta <- c("IAPP", "INS", "DLK1", "INS-IGF2", "G6PC2", 
          "HADH", "ADCYAP1", "GSN", "NPTX2", "C12orf75")

cycling <- c("UBE2C", "TOP2A", "CDK1", "BIRC5", "PBK", "CDKN3",
             "MKI67", "CDC20", "CCNB2", "CDCA3")

immune <- c("ACP5", "APOE", "HLA-DRA", "TYROBP", "LAPTM5", "SDS", 
            "FCER1G", "C1QC", "C1QB")

quiescent_stellate <- c("RGS5", "C11orf96", "FABP4", "CSRP2", "IL24",
                        "ADIRF", "NDUFA4L2", "GPX3", "IGFBP4", "ESAM")

endothelial <- c("PLVAP", "RGCC", "ENG", "PECAM1", "ESM1", "SERPINE1",
                 "CLDN5", "STC1", "MMP1", "GNG11")

schwann <- c("NGFR", "CDH19", "UCN2", "SOX10", "S100A1", "PLP1",
             "TSPAN11", "WNT16", "SOX2", "TFAP2A")

activated_stellate <- c("COL1A1", "COL1A2", "COL6A3", "COL3A1", "TIMP3",
                        "TIMP1", "CTHRC1", "SFRP2", "BGN", "LUM")

epsilon <-c("BHMT", "VSTM2L", "PHGR1", "TM4SF5", "ANXA13", "ASGR1",
            "DEFB1", "GHRL", "COL22A1", "OLFML3")

gamma <- c("PPY", "AQP3", "MEIS2", "ID2", "GPC5-AS1", "CARTPT", "PRSS23",
           "ETV1", "TUBB2A")

delta <- c("SST", "RBP4", "SERPINA1", "RGS2", "PCSK1", "SEC11C", 
           "HHEX", "LEPR", "MDK", "LY6H")

ductal <- c("SPP1", "MMP7", "IGFBP7", "KRT7", "ANXA4", "SERPINA1", "LCN2",
            "CFTR", "KRT19", "SERPING1")

acinar <- c("REG1A", "PRSS1", "CTRB2", "CTRB1", "REG1B", "CELA3A", "PRSS2",
            "REG3A", "CPA1", "CLPS")

alpha <- c("GCG", "TTR", "PPP1R1A", "CRYBA2", "TM4SF4", "MAFB", "GC",
           "GPX3", "PCSK2", "PEMT")

genes <- list(beta, alpha, gamma, delta, 
              acinar, cycling, ductal, immune, 
              activated_stellate, quiescent_stellate, endothelial, epsilon)

azi_markers <- list("Beta" = c("IAPP", "INS", "DLK1", "INS-IGF2", "G6PC2", "HADH", "ADCYAP1", "GSN", "NPTX2", "C12orf75"),
                    "Alpha" = c("GCG", "TTR", "PPP1R1A", "CRYBA2", "TM4SF4", "MAFB", "GC", "GPX3", "PCSK2", "PEMT"),
                    "Delta" = c("SST", "RBP4", "SERPINA1", "RGS2", "PCSK1", "SEC11C", "HHEX", "LEPR", "MDK", "LY6H"),
                    "Gamma" = c("PPY", "AQP3", "MEIS2", "ID2", "GPC5-AS1", "CARTPT", "PRSS23", "ETV1", "TUBB2A"),
                    "Cycling" = c("UBE2C", "TOP2A", "CDK1", "BIRC5", "PBK", "CDKN3", "MKI67", "CDC20", "CCNB2", "CDCA3"),
                    "Acinar" = c("REG1A", "PRSS1", "CTRB2", "CTRB1", "REG1B", "CELA3A", "PRSS2", "REG3A", "CPA1", "CLPS"),
                    "Endothelial" = c("PLVAP", "RGCC", "ENG", "PECAM1", "ESM1", "SERPINE1", "CLDN5", "STC1", "MMP1", "GNG11"),
                    "Q_stellate" = c("RGS5", "C11orf96", "FABP4", "CSRP2", "IL24", "ADIRF", "NDUFA4L2", "IGFBP4", "ESAM"),
                    "A_stellate" = c("COL1A1", "COL1A2", "COL6A3", "COL3A1", "TIMP3", "TIMP1", "CTHRC1", "SFRP2", "BGN", "LUM"),
                    "Immune" = c("ACP5", "APOE", "HLA-DRA", "TYROBP", "LAPTM5", "SDS", "FCER1G", "C1QC", "C1QB", "SRGN"),
                    "Ductal" = c("SPP1", "MMP7", "IGFBP7", "KRT7", "ANXA4", "LCN2", "CFTR", "KRT19", "SERPING1"),
                    "Schwann" = c("NGFR", "CDH19", "UCN2", "SOX10", "S100A1", "PLP1", "TSPAN11", "WNT16", "SOX2", "TFAP2A"))

azi_markers_short <- list("Beta" = c("IAPP", "INS", "NPTX2"),
                    "Alpha" = c("GCG", "TTR", "GC"),
                    "Delta" = c("SST", "RBP4", "SEC11C"),
                    "Gamma" = c("PPY", "MEIS2", "ETV1"),
                    "Cycling" = c("TOP2A", "PBK", "CDKN3"),
                    "Acinar" = c("REG1A", "PRSS1", "CTRB2"),
                    "Endothelial" = c("PLVAP", "RGCC", "ENG"),
                    "Q_stellate" = c("RGS5", "FABP4", "ADIRF"),
                    "A_stellate" = c("COL1A1", "COL1A2", "BGN"),
                    "Immune" = c("ACP5", "APOE", "HLA-DRA"),
                    "Ductal" = c("SPP1", "MMP7", "KRT7"),
                    "Schwann" = c("NGFR", "CDH19", "UCN2"))


purrr::reduce(genes, duplicated)

# Combine all elements into a single vector 
all_elements <- unlist(azi_markers)
# Find duplicated elements 
duplicated_elements <- all_elements[duplicated(all_elements)]
# Remove duplicates from the result 
unique_duplicated_elements <- unique(duplicated_elements)
# Print the duplicated elements 
unique_duplicated_elements

