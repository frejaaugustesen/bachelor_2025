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
motakis_beta@meta.data %>%
  group_by(disease) %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  ungroup() %>%
  ggplot(aes(x = orig.ident, fill = subtype)) +
  geom_bar(position = "fill") +
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

### motakis ----
motakis_beta@meta.data %>%
  group_by(disease) %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  ungroup() %>%
  ggplot(aes(x = disease, fill = subtype)) +
  geom_bar(position = "fill") +
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
