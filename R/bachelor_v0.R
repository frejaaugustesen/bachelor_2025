# Description -------------------------------------------------------------

# Bachelor project spring 2025 at the MadLab group

# first try

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))

set.seed(100)

# installing packages ----------------------------------------------------

## BiocManager 
install.packages("BiocManager")
library(BiocManager) # added to package_load

## scDblFinder - skal bruge hjælp
BiocManager::install("scDblFinder")
library(scDblFinder) # virker ikke??

## Seurat
install.packages("Seurat")
library(Seurat) # added to package_load

## Jointly - skal bruge hjælp
install.packages("remotes")
remotes::install_github("madsen-lab/rJOINTLY")
library(rJOINTLY) # virker ikke??

## qs
install.packages("qs")
library(qs) # added to package_load

## XGBoost - understående er copy paste fra geeksforgeeks.org
# (XGBoost) Installing Packages 
install.packages("data.table") 
install.packages("dplyr") 
install.packages("ggplot2") 
install.packages("caret") 
install.packages("xgboost") 
install.packages("e1071") 
install.packages("cowplot") 

# (XGBoost) Loading packages 
library(data.table) # for reading and manipulation of data 
library(dplyr)     # for data manipulation and joining 
library(ggplot2) # for plotting 
library(caret)     # for modeling 
library(xgboost) # for building XGBoost model 
library(e1071)     # for skewness 
library(cowplot) # for combining multiple plots 

# loading sander seurat object -----------------------------------------------
sander <- qread("data/seurat_objects/wang_sander_merged.qs")


