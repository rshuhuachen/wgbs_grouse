
win_stat <- function(win_start, win_size, funs, values, x){
  
  # subset window data
  x_win <- x %>% filter(pos > win_start & pos <= (win_start+win_size)) %>%
    as.data.frame()
  
  # Vector to hold window stats
  win_stats <- c()
  
  # Calculate statistics for each "values" column
  for(i_col in values){
    n <- sum(!is.na(x_win[, i_col]))
    
    custom_stats <- sapply(funs, function(f, x) eval(parse(text = f))(x), x = x_win[, i_col])
    
    # Return error if functions did not return a single value
    if(is.list(custom_stats)) stop("custom statistics did not return single value")
    
    win_stats <- c(win_stats, n, custom_stats)
  }
  
  return(win_stats)
  
}

win_slide <- function(x, values, win_size, win_step, overlap, funs, cores, region) 
{
  
  # # Warn if window size is bigger than position
  # if(win_size > max(x$pos)){
  #   warning("Window size is bigger than the maximum number of observations in some groups")
  # } 
  # 
  ### Define start position for each window ###
  overlap <- as.numeric(overlap)
  ## if TSS, then only one window (just one value), if no TSS, manual start, middle and end as not equal
  if(region == "TSS"){
    win_start <- min(x$start)
    }  else if (region != "TSS"){
    win_start_first <- c(min(x$start), min(x$start)+overlap)
    win_start_middle <- seq((min(x$start)+overlap+win_size), max(x$end), win_step)
    win_start_middle <- win_start_middle[-(length(win_start_middle))] #minus last one
    win_start_last <- c((max(x$end)-win_size)) 
    win_start <- c(win_start_first, win_start_middle, win_start_last)
  }
  
  ### Loop through windows ###
  funs_out <- parallel::mclapply(win_start, win_stat, win_size, funs, values, x, mc.cores = cores)
  funs_out <- do.call(rbind, funs_out)
  
  ### Make output data.frame ###
  out <- data.frame(win_start = win_start, win_end = win_start + win_size) #took out -1 from original
  out$win_mid <- floor((out$win_start + out$win_end)/2)
  out <- cbind(out, funs_out)
  
  names(out)[4:ncol(out)] <- as.vector(t(outer(values, c("n", funs), paste, sep="_")))
  
  return(out)
}

win_scan <- function(x, 
                    position = NULL, 
                    values, 
                    win_size, 
                    win_step = 0.5*win_size, 
                    overlap,
                    funs = "mean",
                    cores = 1,
                    region){
  ### Parse requested functions to a list ###
  funs <- as.list(funs)
  names(funs) <- sapply(funs, paste)  # name functions' list elements
  
  # ### Check group and position variables ###
  # if(is.null(groups)){
  #   # Create mock group variable
  #   x$group <- 1
  #   groups <- "group"
  # }
  
  if(is.null(position)){
    # Create mock position variable
    x <- x %>% group_by(.dots = groups) %>% mutate(pos = 1:n())
  } else {
    # Rename position variable for downstream functions
    x <- rename(x, "pos" = position)
    assertthat::assert_that(is.numeric(x$pos))
  }
  
  ## Warn if window size is smaller than step size ##
  if(win_size < win_step) warning("Window size is smaller than step size!")
  
  
  ### Compute window statistics per group ###
  out = x %>% 
    do(win_slide(x, values, win_size, win_step, overlap, funs, cores, region)) %>%
    as.data.frame()
  
  return(out)
  
}
