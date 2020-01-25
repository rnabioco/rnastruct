
Sys.setenv(VARNA = "/Users/kriemo/Projects/rnastruct/bin/VARNAv3-93.jar")

#' Plot structure with reactivities (or other values) overlayed using VARNA
#' Reactivites are binned into 10 color bins and colored. Sites without
#' reactivites are not colored. 
#' @param bin If FALSE NA values (and special value -999) will be set to 0 and
#' other values colored using a colorMap. If TRUE, NA values will not be colored, 
#' other values will be colored using a simple color map of 9 red values (white low, red high)
plot_secondary_structure <- function(seq,
                                     ss,
                                     color_by = NULL,
                                     outpath = "tmp.jpeg",
                                     varna_jar = NULL,
                                     bin = TRUE,
                                     windsorize = TRUE,
                                     max_quantile = 0.97,
                                     max_reactivity = NULL,
                                     legend = TRUE,
                                     overlay_args = list(img_pos = c(0.75,
                                                                     0.75))){
  
  if( is.null(varna_jar) && Sys.getenv("VARNA") == "") {
    stop("could not find VARNA environment variable
         please set VARNA=PATH/to/VARNA.jar
         using Sys.setenv(VARNA = 'PATH/to/VARNA.jar'")
  } else {
    varna_jar = Sys.getenv("VARNA") 
  }
  
  structure_str <- ss
  seq_str <- str_replace_all(seq, "T", "U")
  
  output_str <- outpath
  
  varna_args <- c(
    "-cp", varna_jar, "fr.orsay.lri.varna.applications.VARNAcmd",
    "-sequenceDBN", seq_str,
    "-structureDBN", shQuote(structure_str),
    "-zoom", "1", 
    "-resolution", "10.0",
    "-spaceBetweenBases", "0.6",
    "-algorithm", "naview",
    "-o",  output_str
  )
  
  if(!is.null(color_by) && bin){
    # Reds
    col_map <- c("#fff5f0",
                 "#fee0d2",
                 "#fcbba1",
                 "#fc9272",
                 "#fb6a4a",
                 "#ef3b2c",
                 "#cb181d",
                 "#a50f15",
                 "#67000d")
    
    # set white as first color for NA values
    col_map <- c("#ffffff", col_map)
    
    color_by <- as.numeric(color_by)
    color_by[is.na(color_by)] <- -999
    
    reactive_values <- color_by[color_by != -999]
    
    if(length(reactive_values) > 0) {
      if(windsorize){
        if(is.null(max_reactivity)){
          max_value <- quantile(reactive_values, max_quantile)
        } else {
          max_value <- max_reactivity
        }
        reactive_values[reactive_values > max_value] <- max_value
      }
      
      color_max <- max(reactive_values)
      color_min <- min(reactive_values)
    
      breaks <- seq(color_min, color_max, length.out = 9)
      col_breaks <- c(-999.0, breaks)
    } else {
      col_breaks <- c(-999.0)
    }
    
    n_styles <- length(col_map)
    style_args <- paste0("-basesStyle", 1:n_styles, " fill=", col_map, 
                         ",outline=",
                         c("#ffffff", rep("#000000", n_styles - 1)))
    
    base_styles <- findInterval(color_by, col_breaks)
    
    base_style_args <- map_chr(sort(unique(base_styles)),
                               function(x){
                                 paste0("-applyBasesStyle", 
                                        x, 
                                        "on ", 
                                        "'", 
                                        paste0(which(base_styles == x), 
                                               collapse = ","), 
                                        "'")                       
                               })

    varna_args <- c(varna_args, style_args, base_style_args)
    
  } else if (!is.null(color_by)) {
    # RdYlBu
    col_map <- rev(c("#d73027",
                     "#f46d43",
                     "#fdae61",
                     "#fee090",
                     "#ffffbf",
                     "#e0f3f8",
                     "#abd9e9",
                     "#74add1",
                     "#4575b4"))
    
    color_by[color_by == "-999"] <- "0"
    color_str <- paste0(color_by, collapse = ";")
    color_max <- max(as.numeric(color_by))
    color_min <- min(as.numeric(color_by))
    col_breaks <- c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 1) * color_max
    col_map_str <- paste0(col_breaks[-1], ":",col_map[-1], collapse = ";")
    
    colormap_args <- c(
      "-colorMap", shQuote(color_str),
      "-colorMapMax", color_max,
      "-colorMapStyle", shQuote(col_map_str),
      "-colorMapMin", color_min
    )
    
    varna_args <- c(varna_args, colormap_args)
  }
  
  system2("java", 
          varna_args,
          wait = TRUE)
  
  if(legend && bin && !is.null(color_by)){
     f <- tempfile(fileext = ".png")
     on.exit(unlink(f))
    #f <- "tmp.png"
    col_breaks <- col_breaks[col_breaks != -999] 
    col_map <- col_map[col_breaks != -999]
    plot_legend(col_breaks, col_map, fn = f)

    do.call(overlay_images,
            c(png1 = outpath, 
              png2 = f,
              out = outpath,
              overlay_args))
  }
  
}

#' Uses imageMagick to overlay two pngs 
#' Used for adding a legend to a secondary structure
overlay_images <- function(png1, 
                           png2,
                           out = "out.png",
                           img_pos = c(0.5, 0.9)){
  convert_cmd <- Sys.which("convert")
  img_size <- dim(png::readPNG(png1))
  
  h <- img_size[1]
  w <- img_size[2]
  c_args <- c(
    "-size",
    paste0(w,"x",h),
    "-composite",
    png1,
    png2,
    "-geometry",
    paste0(w * 0.5,
           "x",
           h * 0.5,
           "+",
           w * img_pos[1],
           "+",
           h * img_pos[2]),
    "-depth",
    "8",
    out
  )
  
  system2(convert_cmd,
          c_args)
}

#' Generates a legend with continous color map
#' ands saves as png
plot_legend <- function(values, 
                        legend_colors, 
                        fn = "legend.png",
                        ...){
  df <- tibble(x = values, y = values)
  
  p <- ggplot(df, aes(x, y)) +
    geom_point(aes(color = y)) +
    scale_color_gradientn(name = "",
                          colors = legend_colors)

  p <- ggpubr::get_legend(p) %>% ggpubr::as_ggplot()
  
  save_plot(fn, p, base_asp = 1, ...)
}

comp_ss <- function(x, y){
  if(!y %in% c("(", ")")){
    return(0L)
  } else if (x == y){
    return(1L)
  } else {
    return(0L)
  }
}

vcomp_ss <- Vectorize(comp_ss)

calc_ss_similarity <- function(pred, ref, sequence, count_only_ac = TRUE){
  
  lens <- c(nchar(ref), nchar(pred), nchar(sequence))
  if(length(unique(lens)) != 1){
    stop("not all same length")
  }
  
  nts <- str_split(toupper(sequence), "")[[1]]
  
  if(count_only_ac){
    ac_pos = str_which(nts, "A|C")
    sequence <- paste0(nts[ac_pos], collapse = "")
    pred <- paste0(str_split(pred, "")[[1]][ac_pos], collapse = "")
    ref <- paste0(str_split(ref, "")[[1]][ac_pos], collapse = "")
  }
  
  n_bp_ref = str_count(ref, "\\(|\\)")
  n_bp_pred = str_count(pred, "\\(|\\)")
  n_matched_bp = map2_int(str_split(pred, ""),
                          str_split(ref, ""),
                          function(x, y)
                            sum(vcomp_ss(x, y)))
  sen = n_matched_bp / n_bp_ref
  ppv = n_matched_bp / n_bp_pred
  
  sen[is.na(sen)] <- 0
  ppv[is.na(ppv)] <- 0
  
  list(sen = sen,
       ppv = ppv)
}
        


