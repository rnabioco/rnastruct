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

##### Functions ####

#' When writing out excel workbooks using openxlsx::write.xlsx()
#' this function will set the class attributes for a column, which
#' enforces a column type in the resulting xlsx file. 
#' Useful for avoid gene names being clobbered to dates and 
#' setting scientific number formatting

set_xlsx_class <- function(df, col, xlsx_class){
  for(i in seq_along(col)){
    class(df[[col[i]]]) <- xlsx_class
  }
  df
}
