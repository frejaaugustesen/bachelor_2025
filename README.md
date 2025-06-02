# Identification of β-Cell Subtypes in Healthy Individuals and Type 2 Diabetics, Using Machine Learning to Predict Disease States

β-cells are central to glucose homeostasis, and their dysfunction contributes to the pathogenesis of
Type 2 Diabetes (T2D). Recent studies have proposed the existence of distinct β-cell subtypes, with
their proportions potentially shifting across metabolic states. This study aimed to assess the repro-
ducibility and generalizability of β-cell subtypes previously identified by Wang et al (2023), by reim-
plementing their machine learning-based classification approach using single-cell RNA sequencing
(scRNA-seq) data. A classifier was trained on scRNA-seq data from six non-diabetic (ND) and six
T2D donors, as published by Wang et al. (2023), to predict the disease state of individual β-cells based
on gene expression profiles. The model was subsequently applied to β-cells from eight pre-diabetic
(pre-T2D) donors within the same dataset and to an independent scRNA-seq dataset published by
Bandesh et al. (2024), referred to here as the Motakis dataset. Following quality control, 12 donors
from the Motakis dataset were included in the analysis. Both datasets identified two β-cell subtypes,
designated beta1 and beta2. A progressive increase in the proportion of beta2 cells was observed
across ND, pre-T2D, and T2D donors. Additionally, 52 differentially expressed genes distinguishing
the predicted subtypes in the Wang dataset were identified and visualized using heatmaps. Although
the findings suggest transcriptional heterogeneity among β-cells and a possible link to disease pro-
gression, the evidence for distinct, reproducible subtypes remains inconclusive. Further validation
across larger and more diverse cohorts is required to determine these subtypes' robustness and bio-
logical relevance.

![illustration of work](https://github.com/frejaaugustesen/bachelor_2025/blob/main/bachelor_workflow.png)
