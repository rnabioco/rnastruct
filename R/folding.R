
run_rnastructure <- function(seq,
                             ..., 
                             fold_exe = "/Users/kriemo/Projects/rnastruct/bin/RNAstructure/exe/Fold",
                             data_path = "/Users/kriemo/Projects/rnastruct/bin/RNAstructure/data_tables"){
  
  if(Sys.which(fold_exe) == ""){
    stop("RNAStructure Fold not in path")
  }
  
  Sys.setenv(DATAPATH=data_path)
  cmd <- fold_exe
  
  fold_args <- unlist(list(...))
  
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(seq, f)
  
  fold_res <- system2(
    command = cmd,
    args = c(f, "-", "-k", setdiff(fold_args, "-k")),
    stdout = TRUE
  )
  
  fold_res
}


run_duplex_fold <- function(seq1,
                       seq2,
                       ..., 
                       return_dot = TRUE, 
                       fold_exe = "/Users/kriemo/Projects/rnastruct/bin/RNAstructure/exe/DuplexFold",
                       data_path = "/Users/kriemo/Projects/rnastruct/bin/RNAstructure/data_tables"){
  
  if(Sys.which(fold_exe) == ""){
    stop("RNAStructure DuplexFold not in path")
  }
  
  Sys.setenv(DATAPATH=data_path)
  cmd <- fold_exe
  
  fold_args <- unlist(list(...))
  
  f1 <- tempfile()
  on.exit(unlink(f1))
  writeLines(seq1, f1)
  
  
  f2 <- tempfile()
  on.exit(unlink(f2))
  writeLines(seq2, f2)
  
  f3 <- tempfile()
  on.exit(unlink(f3))
  
  fold_res <- system2(
    command = cmd,
    args = c(f1, f2, f3, "-m 1", fold_args),
    stdout = FALSE
  )
  
  if(return_dot){
    res <- ct2dot(f3)
  } else {
    res <- read_delim(f3, trim_ws = TRUE, 
                      delim = " ", 
                      skip = 1, 
                      col_names = c("index",
                                    "nucleotide",
                                    "index1m",
                                    "index1p",
                                    "bp_nt",
                                    "other_idx"),
                      col_types = cols(
                        index = col_double(),
                        nucleotide = col_character(),
                        index1m = col_double(),
                        index1p = col_double(),
                        bp_nt = col_double(),
                        other_idx = col_double()
                      ))
    delG <- read_lines(f3, n_max = 1) %>% 
      str_trim("left") %>% 
      str_split(" +", simplify = TRUE) %>%
      .[4] %>% 
      as.numeric()
    
    res <- list(df = res,
         delg = delG)
  }
  res
}

run_bifold <- function(seq1,
                       seq2,
                       ..., 
                       return_dot = TRUE, 
                       fold_exe = "/Users/kriemo/Projects/rnastruct/bin/RNAstructure/exe/bifold",
                       data_path = "/Users/kriemo/Projects/rnastruct/bin/RNAstructure/data_tables"){
  
  if(Sys.which(fold_exe) == ""){
    stop("RNAStructure bifold not in path")
  }
  
  Sys.setenv(DATAPATH=data_path)
  cmd <- fold_exe
  
  fold_args <- unlist(list(...))
  
  f1 <- tempfile()
  on.exit(unlink(f1))
  writeLines(seq1, f1)
  
  
  f2 <- tempfile()
  on.exit(unlink(f2))
  writeLines(seq2, f2)
  
  f3 <- tempfile()
  on.exit(unlink(f3))
  
  fold_res <- system2(
    command = cmd,
    args = c(f1, f2, f3, "-m 1", fold_args),
    stdout = FALSE
  )
  
  if(return_dot){
    res <- ct2dot(f3)
  } else {
    res <-  tryCatch({
      
    res <- read_delim(f3, trim_ws = TRUE, 
                      delim = " ", 
                      skip = 1, 
                      col_names = c("index",
                                    "nucleotide",
                                    "index1m",
                                    "index1p",
                                    "bp_nt",
                                    "other_idx"),
                      col_types = cols(
                        index = col_double(),
                        nucleotide = col_character(),
                        index1m = col_double(),
                        index1p = col_double(),
                        bp_nt = col_double(),
                        other_idx = col_double()
                      ))
    delG <- read_lines(f3, n_max = 1) %>% 
      str_trim("left") %>% 
      str_split(" +", simplify = TRUE) %>%
      .[4] %>% 
      as.numeric()
    
    res <- list(df = res,
                delg = delG)
    }, error = function(e){
      out <- tribble(
        ~index, ~nucleotide, ~index1m, ~index1p, ~bp_nt, ~other_idx
      )
      res <- list(df = out,
                  delg = NA)
      return(res)
    })
  }
  res
}

run_access_fold <- function(seq1,
                       seq2,
                       ..., 
                       return_dot = TRUE, 
                       fold_exe = "/Users/kriemo/Projects/rnastruct/bin/RNAstructure/exe/AccessFold",
                       data_path = "/Users/kriemo/Projects/rnastruct/bin/RNAstructure/data_tables"){
  
  if(Sys.which(fold_exe) == ""){
    stop("RNAStructure bifold not in path")
  }
  
  Sys.setenv(DATAPATH=data_path)
  cmd <- fold_exe
  
  fold_args <- unlist(list(...))
  
  f1 <- tempfile()
  on.exit(unlink(f1))
  writeLines(seq1, f1)
  
  
  f2 <- tempfile()
  on.exit(unlink(f2))
  writeLines(seq2, f2)
  
  f3 <- tempfile()
  on.exit(unlink(f3))
  
  fold_res <- system2(
    command = cmd,
    args = c(f1, f2, f3, "-m 1", fold_args),
    stdout = FALSE
  )
  
  if(return_dot){
    res <- ct2dot(f3)
  } else {
    res <- read_delim(f3, trim_ws = TRUE, 
                      delim = " ", 
                      skip = 1, 
                      col_names = c("index",
                                    "nucleotide",
                                    "index1m",
                                    "index1p",
                                    "bp_nt",
                                    "other_idx"),
                      col_types = cols(
                        index = col_double(),
                        nucleotide = col_character(),
                        index1m = col_double(),
                        index1p = col_double(),
                        bp_nt = col_double(),
                        other_idx = col_double()
                      ))
  }
  res
}

ct2dot <- function(fn, 
                   ct2dot_exe = "/Users/kriemo/Projects/rnastruct/bin/RNAstructure/exe/ct2dot"){
  if(Sys.which(ct2dot_exe) == ""){
    stop("RNAStructure ct2dot not in path")
  }
  res <- system2(ct2dot_exe,
          args = c(fn, "1", "-"),
          stdout = TRUE
          )
  res
}

nts <- c("A", "T", "C", "G")
randomize_seq <- function(n){
  sample(nts, n, replace = TRUE)
}

shuffle_seq <- function(seq, seed = NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }
  n_bps <- nchar(seq)
  rand_idx <- sample(1:n_bps, n_bps, replace = FALSE)
  rand_seq <- str_split(seq, "", simplify = T)[rand_idx]
  str_c(rand_seq, collapse = "")
}

shuffle_seqs <- function(seqs){
  map_chr(seqs, shuffle_seq)
}


### 











