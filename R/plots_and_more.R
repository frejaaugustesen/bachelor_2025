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


# save --------------------------------------------------------------------
qsave(motakis_beta, file = here::here("data/seurat_objects/motakis_beta_new.qs"))



# renaming subtypes -------------------------------------------------------


wang_beta@meta.data <- wang_beta@meta.data %>%
  mutate(subtype = dplyr::recode(subtype,
                          "nd" = "beta1",
                          "t2d" = "beta2"))


motakis_beta@meta.data <- motakis_beta@meta.data %>%
  mutate(subtype = dplyr::recode(subtype,
                          "nd" = "beta1",
                          "t2d" = "beta2"))

# Stacked barplots --------------------------------------------------------

# laver om til procent
motakis_tabel <- motakis_beta@meta.data %>% 
  as.data.frame() %>% 
  group_by(disease, subtype, orig.ident, sample) %>% 
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

# laver om til procent
wang_tabel <- wang_beta@meta.data %>% 
  as.data.frame() %>% 
  group_by(disease, subtype, orig.ident) %>% 
  tally() %>% 
  group_by(orig.ident) %>% 
  mutate(perc = (n/sum(n)*100), 
         sum = sum(n))

wang_tabel_disease <- wang_beta@meta.data %>% 
  as.data.frame() %>% 
  group_by(disease, subtype) %>% 
  tally() %>% 
  group_by(disease) %>% 
  mutate(perc = (n/sum(n)*100), 
         sum = sum(n))

table(wang_beta@meta.data$subtype, wang_beta@meta.data$disease)

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
    "beta1" = "lightblue",       
    "beta2" = "pink"     
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
    "beta1" = "lightblue",       
    "beta2" = "pink"     
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
    title = "Predicted subtypes per donor in Wang Sander",
    x = "Disease",
    y = "Percentage"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.y = element_line(color = "grey"),  
        axis.ticks.y = element_line(color = "grey"),
        axis.ticks.x = element_blank())+
  scale_fill_manual(values = c(
    "beta1" = "lightblue",       
    "beta2" = "pink"     
  ))

### Motakis ----
motakis_beta@meta.data %>%
  group_by(disease) %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  ungroup() %>%
  ggplot(aes(x = disease, fill = subtype)) +
  geom_bar(position = "fill") +
  labs(
    title = "Predicted subtypes per donor in Motakis",
    x = "Disease",
    y = "Percentage"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.y = element_line(color = "grey"),  
        axis.ticks.y = element_line(color = "grey"),
        axis.ticks.x = element_blank())+
  scale_fill_manual(values = c(
    "beta1" = "lightblue",       
    "beta2" = "pink"     
  ))


# 2c like figures ---------------------------------------------------------

## wang ----

wang_beta1 <- wang_beta1_donor %>% select(disease, perc) %>% mutate(subtype = "beta1")
wang_beta2 <- wang_beta2_donor %>% select(disease, perc) %>% mutate(subtype = "beta2")

wang_combined_df <- bind_rows(wang_beta1, wang_beta2)

wang_tabel_disease %>% 
  ggplot(aes(x = subtype, y = perc, fill = disease)) +  
  geom_col(position = position_dodge(width = 0.9)) +
  ylim(0, 100) +
  geom_point(data = wang_combined_df, 
             aes(x= subtype, y = perc, group = disease), 
             position = position_dodge(width = 0.9),
             size = 1, color = "darkgrey", alpha = 0.6) +
  labs(
    title = "Subtype distribution in Wang Sander",
    x = "Subtype",
    y = "Percentage",
    fill = "Disease"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    text = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.y = element_line(color = "grey"),  
    axis.ticks.y = element_line(color = "grey")
  ) +
  scale_fill_manual(values = c(
    "nd" = "lightblue",       
    "pre" = "lightgreen",      
    "t2d" = "pink"     
  )) +
  # Add manual p-value labels
  annotate("text", x = 0.85, y = 80, label = "P = 0.002", size = 4) +
  annotate("text", x = 0.85, y = 76, label = "**", size = 5) +
  annotate("text", x = 1.15, y = 45, label = "P = 0.061", size = 4) +
  annotate("text", x = 1.85, y = 75, label = "P = 0.002", size = 4) +
  annotate("text", x = 1.85, y = 71, label = "**", size = 5) +
  annotate("text", x = 2.15, y = 95, label = "P = 0.061", size = 4)

## motakis ----

motakis_beta1 <- motakis_tabel %>% select(disease, perc) %>% mutate(subtype = "beta1")
motakis_beta2 <- motakis_tabel %>% select(disease, perc) %>% mutate(subtype = "beta2")


motakis_beta1 <- motakis_beta1_donor %>% select(disease, perc) %>% mutate(subtype = "beta1")
motakis_beta2 <- motakis_beta2_donor %>% select(disease, perc) %>% mutate(subtype = "beta2")

motakis_combined_df <- bind_rows(motakis_beta1, motakis_beta2)

motakis_tabel_disease %>% 
  ggplot(aes(x = subtype, y = perc, fill = disease)) +  
  geom_col(position = position_dodge(width = 0.9)) +
  ylim(0, 100) +
  geom_point(data = motakis_combined_df, 
             aes(x= subtype, y = perc, group = disease), 
             position = position_dodge(width = 0.9),
             size = 1, color = "darkgrey", alpha = 0.6) +
  labs(
    title = "Subtype distribution in Motakis",
    x = "Subtype",
    y = "Percentage",
    fill = "Disease"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    text = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.y = element_line(color = "grey"),  
    axis.ticks.y = element_line(color = "grey")
  ) +
  scale_fill_manual(values = c(
    "nd" = "lightblue",       
    "pre" = "lightgreen",      
    "t2d" = "pink"     
  )) +
  # Add manual p-value labels
  annotate("text", x = 0.85, y = 72, label = "P = 0.954", size = 4) +
  annotate("text", x = 1.15, y = 63, label = "P = 0.476", size = 4) +
  annotate("text", x = 1.85, y = 51, label = "P = 0.954", size = 4) +
  annotate("text", x = 2.15, y = 60, label = "P = 0.476", size = 4)


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

