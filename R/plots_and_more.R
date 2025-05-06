# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group
# making plots for the thesis

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))

set.seed(100)


# load --------------------------------------------------------------------

motakis_beta <- qread("/work/bachelor_2025/data/seurat_objects/motakis_beta_integrated.qs")
wang_beta <- qread("/work/bachelor_2025/data/seurat_objects/wang_beta_pred.qs")


# Stacked barplots --------------------------------------------------------

# laver om til procent
motakis_tabel <- motakis_beta@meta.data %>% 
  as.data.frame() %>% 
  group_by(disease, subtype, orig.ident) %>% 
  tally() %>% 
  group_by(orig.ident) %>% 
  mutate(perc = (n/sum(n)*100), 
         sum = sum(n))

motakis_tabel_disease <- motakis_beta@meta.data %>% 
  as.data.frame() %>% 
  group_by(disease, subtype) %>% 
  tally() %>% 
  group_by(disease) %>% 
  mutate(perc = (n/sum(n)*100), 
         sum = sum(n))

table(motakis_beta@meta.data$subtype, motakis_beta@meta.data$disease)

## per donor --------------------------------------------------------------

### wang ----
wang_beta@meta.data %>%
  group_by(disease) %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  ungroup() %>%
  ggplot(aes(x = orig.ident, fill = subtype)) +
  geom_bar(position = "fill") +
  facet_wrap(~ disease, scales = "free_x") +
  labs(
    title = "Predicted subtypes per donor (Wang Sander)",
    x = "Donor",
    y = "Percentage"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c(
    "nd" = "lightblue",       
    "t2d" = "pink"     
  ))

### motakis ----
motakis_tabel %>%
  group_by(disease) %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  ungroup() %>%
  ggplot(aes(x = orig.ident, y = perc, fill = subtype)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_wrap(~ disease, scales = "free_x") +
  labs(
    title = "Predicted subtypes per donor (Motakis)",
    x = "Donor",
    y = "Percentage"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c(
    "nd" = "lightblue",       
    "t2d" = "pink"     
  ))


## per disease -------------------------------------------------------------

### wang ----
wang_beta@meta.data %>%
  group_by(disease) %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  ungroup() %>%
  ggplot(aes(x = disease, fill = subtype)) +
  geom_bar(position = "fill") +
  labs(
    title = "Predicted subtypes per donor (Wang Sander)",
    x = "Donor",
    y = "Percentage"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c(
    "nd" = "lightblue",       
    "t2d" = "pink"     
  ))

motakis_beta@meta.data %>%
  group_by(disease) %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  ungroup() %>%
  ggplot(aes(x = disease, fill = subtype)) +
  geom_bar(position = "fill") +
  labs(
    title = "Predicted subtypes per donor (Wang Sander)",
    x = "Donor",
    y = "Percentage"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c(
    "nd" = "lightblue",       
    "t2d" = "pink"     
  ))

### motakis ----
motakis_tabel_disease %>%
  #group_by(disease) %>%
  #mutate(orig.ident = factor(orig.ident)) %>%
  ungroup() %>%
  ggplot(aes(x = disease, y=perc, fill = subtype)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_point(data = motakis_tabel, aes(x= disease, y = perc, color = subtype,
                                       alpha = 0.5), 
             position = position_dodge(width = 1)) +
  labs(
    title = "Predicted subtypes per donor (Motakis)",
    x = "Donor",
    y = "Percentage"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(
    "nd" = "lightblue",       
    "t2d" = "pink"     
  ))


# 2c like figures ---------------------------------------------------------

## wang ----
wang_beta@meta.data %>%
  ggplot(aes(x = subtype, fill = disease)) +  
  geom_bar(position = position_dodge(width = 0.9)) +
  labs(
    title = "Relative abundance of beta cell subtype in ND, Pre and T2D (Wang Sander)",
    x = "Subtype",
    y = "Cell count",
    fill = "Disease"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    text = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_fill_manual(values = c(
    "nd" = "lightblue",       
    "pre" = "lightgreen",      
    "t2d" = "pink"     
  ))

## motakis ----
motakis_beta@meta.data %>%
  ggplot(aes(x = subtype, fill = disease)) +  
  geom_bar(position = position_dodge(width = 0.9)) +
    geom_point()+
  labs(
    title = "Relative abundance of beta cell subtype in ND, Pre and T2D (Motakis)",
    x = "Subtype",
    y = "Cell count",
    fill = "Disease"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    text = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_fill_manual(values = c(
    "nd" = "lightblue",       
    "pre" = "lightgreen",      
    "t2d" = "pink"     
  ))


# wang study comparison (UMAP) ---------------------------------------------------

DimPlot(wang_beta, reduction = "umap", label = TRUE,
        group.by = c("subtype", "study_subtype", "agreement"),
        label.size = 2.5,
        repel = TRUE) &
  NoLegend()

DimPlot(wang_beta, reduction = "umap", label = TRUE,
        group.by = "agreement",
        label.size = 2.5,
        repel = TRUE)

