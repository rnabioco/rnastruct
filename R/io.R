

read_dp <- function(fn, ...){
  x <- read_lines(fn, skip = 1, ...)
  tibble(seq = x[1],
         ss = x[2])
}

read_db <- read_dp
#' read .map file into df
read_map <- function(fn){
  x <- read_tsv(fn, 
                col_names = c("idx", "react", "stderr", "nuc"),
                col_types = c("iddc"))
  x
}

read_rnastructure <- function(fn){
  
  lines <- readr::read_lines(fn)
  first_char <- str_sub(lines, 1, 1)
  hdr_idx <- which(first_char == ">")
  ss_idx <- which(first_char %in% c(".", "(", ")"))
  
  df <- tibble(header = lines[hdr_idx],
               ss = lines[ss_idx])
  
  df <- mutate(df,
               header = str_remove(header, "^>"),
               header = str_trim(header),
               header = str_split(df$header, " ", simplify = F) %>%
                 map_chr(~.x[length(.x)]),
               header = str_remove(header, "^>"))
  df
}


read_rnafold <- function(fn) {
  x <- read_lines(fn)
  
  first_char <- str_sub(x, 1, 1)
  hdr_idx <- seq(1, length(x), by = 3)
  seq_idx <- seq(2, length(x), by = 3)
  ss_idx <- seq(3, length(x), by = 3)
  
  res <- tibble(header = x[hdr_idx],
                seq = x[seq_idx],
                ss = x[ss_idx])
  res <- mutate(res, 
                ss = str_split(ss, " ", simplify = TRUE) %>% .[, 1],
                header = str_remove(header, "^>"))
  res
}





read_tabix <- function(tbx, region_df = NULL, region = NULL,
                       file_query = "-R",
                       group_by_query = FALSE){
  
  tabix <-  Sys.which("tabix")
  
  if(tabix == ""){
    stop("couldn't find tabix")
  }
  
  if(!is.null(region_df) && !is.null(region)){
    stop("please supply either region_df or region, not both")
  }
  
  if(!is.null(region_df)){
    cols_for_region_df <- c("chrom", "start", "end")
    
    if(!all(cols_for_region_df %in% colnames(region_df))){
      stop("missing columns in region_df: ", 
           paste0(cols_for_region_df, collapse = " "))
    }
    
    res <- format_tsv(region_df[cols_for_region_df], 
                      col_names = FALSE) %>% 
      system2(tabix, 
              c(tbx, 
                file_query,
                "-"),
              stdout = TRUE,
              input = .) 
  } else if (!is.null(region)) {
    res <- system2(tabix, 
                   c(tbx,
                     region),
                   stdout = TRUE) 
  }
  
  if (is.null(region) && is.null(region_df)){
    res <- gzfile(tbx)
  }
  
  if(length(res) == 0){
    res <- tibble()
  } else {
    res <- read_tsv(res, col_names = FALSE)
  }
  
  res
}



read_bw <- function(bw_fn,
                   chrom = NULL,
                   start = NULL,
                   end = NULL,
                   bigWigToBedGraph = "/Users/kriemo/bin/kent/bigWigToBedGraph") {
  
  ucsc_args <- vector("character")
  
  if (!is.null(chrom)) {
    ucsc_args <- c(ucsc_args, paste0("-chrom=", chrom))
  }
  
  if (!is.null(start)) {
    ucsc_args <- c(ucsc_args, paste0("-start=", start))
  }
  
  if (!is.null(end)) {
    ucsc_args <- c(ucsc_args, paste0("-end=", end))
  }
  
  ucsc_args <- c(ucsc_args,
                 bw_fn,
                 "stdout")
  
  res <- read_tsv(
    system2(bigWigToBedGraph,
            ucsc_args,
            stdout = TRUE),
    col_names = c("chrom",
                  "start",
                  "end",
                  "coverage")
  )
  
  res
}
