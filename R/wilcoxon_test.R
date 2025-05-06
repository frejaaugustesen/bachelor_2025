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


# renaming subtypes -------------------------------------------------------

wang_beta@meta.data <- wang_beta@meta.data %>%
  mutate(subtype = recode(subtype,
                          "nd" = "beta1",
                          "t2d" = "beta2"))

motakis_beta@meta.data <- motakis_beta@meta.data %>%
  mutate(subtype = recode(subtype,
                          "nd" = "beta1",
                          "t2d" = "beta2"))

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


# ANOVA -------------------------------------------------------------------
# performing one-way ANOVA test on percentage of subtypes in disease states

# Null hypothesis: the means of the different groups are the same
# Alternative hypothesis: At least one sample mean is not equal to the others.

## assumptions rules ----
# The observations are obtained independently and randomly from the population 
# defined by the factor levels

# The data of each factor level are normally distributed.

# These normal populations have a common variance. (Levene’s test can be used
# to check this.) see this further down

## Wang ----
wang_tabel_disease <- wang_beta@meta.data %>% 
  as.data.frame() %>% 
  group_by(disease, subtype) %>% 
  tally() %>% 
  group_by(disease) %>% 
  mutate(perc = (n/sum(n)*100), 
         sum = sum(n))

wang_tabel <- wang_beta@meta.data %>% 
  as.data.frame() %>% 
  group_by(disease, subtype, orig.ident) %>% 
  tally() %>% 
  group_by(orig.ident) %>% 
  mutate(perc = (n/sum(n)*100), 
         sum = sum(n))

### structure of data ----
# dividing into to dataframes - one for beta1 and one for beta2

# Subset for beta1
wang_beta1 <- wang_tabel_disease %>%
  filter(subtype == "beta1") %>% 
  as.data.frame()

wang_beta1_donor <- wang_tabel %>%
  filter(subtype == "beta1") %>% 
  as.data.frame()

wang_beta1_donor$disease <- as.factor(wang_beta1_donor$disease)

# Subset for beta2
wang_beta2 <- wang_tabel_disease %>%
  filter(subtype == "beta2") %>% 
  as.data.frame()

wang_beta2_donor <- wang_tabel %>%
  filter(subtype == "beta2") %>% 
  as.data.frame()

wang_beta2_donor$disease <- as.factor(wang_beta1_donor$disease)

# Now levels will work
levels(wang_beta1_donor$disease)
levels(wang_beta2_donor$disease)

# performing assumption about normal populations have a common variance 
# using Levene’s test
leveneTest(perc ~ disease, data = wang_beta1_donor)
leveneTest(perc ~ disease, data = wang_beta2_donor)

### plot ----
# plot to visualize and help with understanding

wang_tabel_disease %>% 
  ggplot(aes(x = subtype, y = perc, fill = disease)) +  
  geom_col(position = position_dodge(width = 0.9)) +
  labs(
    title = "Relative abundance of beta cell subtype in ND, Pre and T2D (Wang Sander)",
    x = "Subtype",
    y = "Percentage",
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

### Check data ----
wang_beta1_donor 
wang_beta2_donor

levels(wang_beta1_donor$disease)
levels(wang_beta2_donor$disease)
# levels are in aplhabetic order, nice!

# Compute summary statistics by groups - count, mean, sd
group_by(wang_beta1_donor, disease) %>%
  summarise(
    count = n(),
    mean = mean(perc, na.rm = TRUE),
    sd = sd(perc, na.rm = TRUE)
  )

group_by(wang_beta2_donor, disease) %>%
  summarise(
    count = n(),
    mean = mean(perc, na.rm = TRUE),
    sd = sd(perc, na.rm = TRUE)
  )

### Visualize data ----

# Box plots
# Plot perc by disease and color by disease

# beta1
ggboxplot(wang_beta1_donor, x = "disease", y = "perc", 
          color = "disease", palette = c("lightblue", "lightgreen", "pink"),
          order = c("nd", "pre", "t2d"),
          ylab = "perc", xlab = "disease")

#beta2
ggboxplot(wang_beta2_donor, x = "disease", y = "perc", 
          color = "disease", palette = c("lightblue", "lightgreen", "pink"),
          order = c("nd", "pre", "t2d"),
          ylab = "perc", xlab = "disease")

# Mean plots
# Plot perc by disease
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)

# beta1
ggline(wang_beta1_donor, x = "disease", y = "perc", 
       add = c("mean_se", "jitter"), 
       order = c("nd", "pre", "t2d"),
       ylab = "perc", xlab = "Disease")

# beta2
ggline(wang_beta2_donor, x = "disease", y = "perc", 
       add = c("mean_se", "jitter"), 
       order = c("nd", "pre", "t2d"),
       ylab = "perc", xlab = "Disease")

### Anova-test ----

# Compute the analysis of variance
wang_res.aov_beta1 <- aov(perc ~ disease, data = wang_beta1_donor)
# Summary of the analysis
summary(wang_res.aov_beta1)

# Compute the analysis of variance
wang_res.aov_beta2 <- aov(perc ~ disease, data = wang_beta2_donor)
# Summary of the analysis
summary(wang_res.aov_beta2)

### Tukey HSD ----
# Tukey Honest Significant Differences
# used to determine if there is significant difference between the specific groups

TukeyHSD(wang_res.aov_beta1)
TukeyHSD(wang_res.aov_beta2)
# not significant difference between pre and t2d!

# other test??
pairwise.t.test(wang_beta1_donor$perc, wang_beta1_donor$disease,
                p.adjust.method = "BH")

### assumptions check ----

# 1. Homogeneity of variances
plot(res.aov_beta1, 1)
# 1. Homogeneity of variances
plot(res.aov_beta2, 1)

# 2. Normality
plot(res.aov_beta1, 2)
plot(res.aov_beta2, 2)

# Extract the residuals
aov_residuals_beta1 <- residuals(object = res.aov_beta1 )
aov_residuals_beta2 <- residuals(object = res.aov_beta2 )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals_beta1 )
shapiro.test(x = aov_residuals_beta2 )

### new plot ----

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
    title = "Relative abundance of beta cell subtype in ND, Pre and T2D (Wang Sander)",
    x = "Subtype",
    y = "Percentage",
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
  )) +
  # Add manual p-value labels
  annotate("text", x = 0.85, y = 80, label = "P = 0.002", size = 3) +
  annotate("text", x = 0.85, y = 76, label = "**", size = 4) +
  annotate("text", x = 1.15, y = 45, label = "P = 0.061", size = 3) +
  annotate("text", x = 1.85, y = 75, label = "P = 0.002", size = 3) +
  annotate("text", x = 1.85, y = 71, label = "**", size = 4) +
  annotate("text", x = 2.15, y = 95, label = "P = 0.061", size = 3)


## Motakis ----
motakis_tabel_disease <- motakis_beta@meta.data %>% 
  as.data.frame() %>% 
  group_by(disease, subtype) %>% 
  tally() %>% 
  group_by(disease) %>% 
  mutate(perc = (n/sum(n)*100), 
         sum = sum(n))

motakis_tabel <- motakis_beta@meta.data %>% 
  as.data.frame() %>% 
  group_by(disease, subtype, orig.ident) %>% 
  tally() %>% 
  group_by(orig.ident) %>% 
  mutate(perc = (n/sum(n)*100), 
         sum = sum(n))

### structure of data ----
# dividing into to dataframes - one for beta1 and one for beta2

# Subset for beta1
motakis_beta1 <- motakis_tabel_disease %>%
  filter(subtype == "beta1") %>% 
  as.data.frame()

motakis_beta1_donor <- motakis_tabel %>%
  filter(subtype == "beta1") %>% 
  as.data.frame()

motakis_beta1_donor$disease <- as.factor(motakis_beta1_donor$disease)

# Subset for beta2
motakis_beta2 <- motakis_tabel_disease %>%
  filter(subtype == "beta2") %>% 
  as.data.frame()

motakis_beta2_donor <- motakis_tabel %>%
  filter(subtype == "beta2") %>% 
  as.data.frame()

motakis_beta2_donor$disease <- as.factor(motakis_beta1_donor$disease)

# Now levels will work
levels(motakis_beta1_donor$disease)
levels(motakis_beta2_donor$disease)

# performing assumption about normal populations have a common variance 
# using Levene’s test
leveneTest(perc ~ disease, data = motakis_beta1_donor)
leveneTest(perc ~ disease, data = motakis_beta2_donor)

### plot ----
# plot to visualize and help with understanding

motakis_tabel_disease %>% 
  ggplot(aes(x = subtype, y = perc, fill = disease)) +  
  geom_col(position = position_dodge(width = 0.9)) +
  labs(
    title = "Relative abundance of beta cell subtype in ND, Pre and T2D (motakis Sander)",
    x = "Subtype",
    y = "Percentage",
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

### Check data ----
motakis_beta1_donor 
motakis_beta2_donor

levels(motakis_beta1_donor$disease)
levels(motakis_beta2_donor$disease)
# levels are in aplhabetic order, nice!

# Compute summary statistics by groups - count, mean, sd
group_by(motakis_beta1_donor, disease) %>%
  summarise(
    count = n(),
    mean = mean(perc, na.rm = TRUE),
    sd = sd(perc, na.rm = TRUE)
  )

group_by(motakis_beta2_donor, disease) %>%
  summarise(
    count = n(),
    mean = mean(perc, na.rm = TRUE),
    sd = sd(perc, na.rm = TRUE)
  )

### Visualize data ----

# Box plots
# Plot perc by disease and color by disease

# beta1
ggboxplot(motakis_beta1_donor, x = "disease", y = "perc", 
          color = "disease", palette = c("lightblue", "lightgreen", "pink"),
          order = c("nd", "pre", "t2d"),
          ylab = "perc", xlab = "disease")

#beta2
ggboxplot(motakis_beta2_donor, x = "disease", y = "perc", 
          color = "disease", palette = c("lightblue", "lightgreen", "pink"),
          order = c("nd", "pre", "t2d"),
          ylab = "perc", xlab = "disease")

# Mean plots
# Plot perc by disease
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)

# beta1
ggline(motakis_beta1_donor, x = "disease", y = "perc", 
       add = c("mean_se", "jitter"), 
       order = c("nd", "pre", "t2d"),
       ylab = "perc", xlab = "Disease")

# beta2
ggline(motakis_beta2_donor, x = "disease", y = "perc", 
       add = c("mean_se", "jitter"), 
       order = c("nd", "pre", "t2d"),
       ylab = "perc", xlab = "Disease")

### Anova-test ----

# Compute the analysis of variance
motakis_res.aov_beta1 <- aov(perc ~ disease, data = motakis_beta1_donor)
# Summary of the analysis
summary(motakis_res.aov_beta1)

# Compute the analysis of variance
motakis_res.aov_beta2 <- aov(perc ~ disease, data = motakis_beta2_donor)
# Summary of the analysis
summary(motakis_res.aov_beta2)

### Tukey HSD ----
# Tukey Honest Significant Differences
# used to determine if there is significant difference between the specific groups

TukeyHSD(motakis_res.aov_beta1)
TukeyHSD(motakis_res.aov_beta2)
# not significant difference between pre and t2d!

# other test??
pairwise.t.test(motakis_beta1_donor$perc, motakis_beta1_donor$disease,
                p.adjust.method = "BH")

### assumptions check ----

# 1. Homogeneity of variances
plot(res.aov_beta1, 1)
# 1. Homogeneity of variances
plot(res.aov_beta2, 1)

# 2. Normality
plot(res.aov_beta1, 2)
plot(res.aov_beta2, 2)

# Extract the residuals
aov_residuals_beta1 <- residuals(object = res.aov_beta1 )
aov_residuals_beta2 <- residuals(object = res.aov_beta2 )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals_beta1 )
shapiro.test(x = aov_residuals_beta2 )

### new plot ----

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
    title = "Relative abundance of beta cell subtype in ND, Pre and T2D (Wang Sander)",
    x = "Subtype",
    y = "Percentage",
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
  )) +
  # Add manual p-value labels
  annotate("text", x = 0.85, y = 68, label = "P = 0.735", size = 3) +
  annotate("text", x = 1.15, y = 58, label = "P = 0.537", size = 3) +
  annotate("text", x = 1.85, y = 58, label = "P = 0.735", size = 3) +
  annotate("text", x = 2.15, y = 68, label = "P = 0.537", size = 3)


