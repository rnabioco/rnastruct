# libs and constants shared across markdown docs

library(valr)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(Matrix)
library(R.utils)
library(viridis)
library(jsonlite)
library(rtracklayer)
library(kentr)
library(ComplexHeatmap)

#### Paths ####

project_dir <- path.expand("~/Projects/rnastruct")
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
docs_dir <- file.path(project_dir, "docs")
db_dir <- file.path(project_dir, "dbases")

theme_set(theme_cowplot())
