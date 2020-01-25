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

#### annotation files ####

annots <- list(
  fa = "~/Projects/shared_dbases/genomes/human/chr_appended_grch37/Homo_sapiens.GRCh37.dna.primary_assembly_chrappended.fa",
  fa_idx = "~/Projects/shared_dbases/genomes/human/chr_appended_grch37/Homo_sapiens.GRCh37.dna.primary_assembly_chrappended.fa.fai"
)

#### Other R scripts ####
source_files <- dir(file.path(project_dir, "R"), 
                    pattern = ".R$",
                    full.names = TRUE)

source_files <- setdiff(source_files, 
                        file.path(project_dir, "R", "constants.R"))

walk(source_files, source)
