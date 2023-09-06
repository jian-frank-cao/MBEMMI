## Functions to clean raw data -------------------------------------------------
# Function to normalize the data
normalize <- function(data) {
  
  # Subtract mean and divide by standard deviation for each column
  normalized_data <- data %>%
    mutate(across(everything(), 
                  ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)
    ))
  
  # Replace NA values with 0                  ### NA may be caused by sd = 0  
  normalized_data[is.na(normalized_data)] <- 0
  
  return(normalized_data)
}


# Prepare auxiliary variables
prepare_auxiliary <- function(auxiliary) {
  
  # Check if auxiliary is not null
  if (!is.null(auxiliary)) {
    
    # Normalize and rename columns
    vars_aux <- auxiliary %>%
      normalize() %>%
      dplyr::rename_with(~paste0("auxiliary-", 1:length(.)))
    
  } else {
    vars_aux <- NULL
  }
  
  return(vars_aux)
}


# Create polytime variables
create_polytime <- function(n_time, n_polytime) {
  
  # Generate polytime matrix and normalize
  vars_polytime <- 1:n_polytime %>%
    purrr::map(~ (1:n_time) ^ .) %>%
    do.call("cbind", .) %>%
    as.data.frame() %>%
    normalize() %>%
    dplyr::rename_with(~ paste0("polytime-", 1:n_polytime))
  
  return(vars_polytime)
}


# Clean data
clean_data <- function(data, n_polytime = 3){
  n_industry <-  ncol(data$detail)
  n_time <-  nrow(data$detail)
  
  # Rescale data
  multiplier <- mean(colMeans(data$detail, na.rm = TRUE), na.rm = TRUE)
  data$multiplier <- multiplier
  data$detail <- data$detail / multiplier
  data$q_total <- data$q_total / multiplier
  data$a_total <- data$a_total / multiplier
  data$a_q_total <- data$a_q_total / multiplier
  
  # Rename columns of 'detail'
  colnames(data$detail) <- paste0("industry-", 1:n_industry)
  
  # Add auxiliary variables
  if (!is.null(data$auxiliary)) {
    data$detail_with_aux <- cbind(
      data$detail,
      prepare_auxiliary(data$auxiliary),
      create_polytime(n_time, n_polytime)
    )
  }else{
    data$detail_with_aux <- cbind(
      data$detail,
      create_polytime(n_polytime)
    )
  }
  
  # Add intercept and reorder columns to put 'intercept' first
  data$detail_with_aux <- data$detail_with_aux %>%
    mutate(intercept = 1) %>%
    dplyr::select(intercept, everything())
  
  return(data)
}

## Functions to prepare Sobol sequence -----------------------------------------
# Prepare Sobol squence
prepare_sobol <- function(dim, n = 1000, scrambling = 0, seed = 316){
  
  # Randomly decide on a 'burn-in' period between 0 and 1000.
  burn_in <- sample(0:1000, 1)
  
  # Generate the Sobol sequence
  sobol_seq <- randtoolbox::sobol(
    n = n + burn_in,
    dim = dim,
    init = TRUE,
    scrambling = scrambling,
    seed = seed
  )
  
  # Skipping the 'burn-in' period
  sobol_seq <- sobol_seq[(burn_in + 1):(n + burn_in), ]
  
  return(sobol_seq)
}


## Functions to do obs-by-obs bootstrapphing -----------------------------------
# bootstrapping = function(n_row, n_bootstrap, 
#                          sobol_seq = NULL, qmc = TRUE,
#                          scrambling = 0, seed = 316){
#   # sobol seq
#   if (is.null(sobol_seq)) {
#     sobol_seq = prepare_sobol(
#       dim = n_bootstrap,
#       n = n_row,
#       scrambling = scrambling,
#       seed = seed
#     )
#   }
#   
#   ind_boots = NULL
#   for (i in 1:n_bootstrap) {
#     if (qmc == TRUE) {
#       # setup
#       qmc_seq = sobol_seq[,i]
#       prob = 1/nrows
#       n_picked = 0
#       finished = FALSE
#       times_picked = rep(0, nrows)
#       
#       # determine # of times being picked
#       for (row in 1:n_row) {
#         if (!finished) {
#           times_could_be_picked = 0:(n_row - n_picked)
#           pbinom_prob = prob/(1-(row-1)*prob)
#           if (pbinom_prob > 1) {
#             pbinom_prob = 1
#           }
#           y = pbinom(
#             q = times_could_be_picked,
#             size = n_row - n_picked,
#             prob = pbinom_prob
#           )
#           times_picked[row] = (y > qmc_seq[row]) %>% 
#             which %>% 
#             min - 1
#           n_picked = n_picked + times_picked[row]
#         }
#         if (n_picked == n_row) {
#           finished = TRUE
#         }
#       }
#       
#       # shuffle the picked rows
#       ind = mapply(rep, 1:nrows, times_picked) %>% 
#         do.call("c", .) %>% 
#         sample(., nrows, replace = FALSE)
#     } else {
#       ind = sample(1:nrows, nrows, replace = TRUE)
#     }
#     
#     # get indices for original to bootstrapped
#     orig_boots = ind
#     
#     # get indices for boostrapped to original
#     boots_orig = match(1:length(ind), ind)
#     
#     ind_boots[[paste0("bootstrap-", i)]] = list(
#       orig_boots = orig_boots,
#       boots_orig = boots_orig
#     )
#   }
#   
#   # finish
#   output = list(
#     ind_boots = ind_boots,
#     sobol_seq = sobol_seq
#   )
#   
#   return(output)
# }


## Functions to perform observation-by-observation bootstrapping
bootstrapping <- function(n_row, n_bootstrap, 
                          sobol_seq = NULL, qmc = TRUE,
                          scrambling = 0, seed = 316){
  
  # If no Sobol sequence is provided, generate one
  if (is.null(sobol_seq)) {
    sobol_seq <- prepare_sobol(
      dim = n_bootstrap,
      n = n_row,
      scrambling = scrambling,
      seed = seed
    )
  }
  
  ind_boots <- list()
  
  for (i in 1:n_bootstrap) {
    
    # If Quasi-Monte Carlo method is selected
    if (qmc) {
      
      qmc_seq <- sobol_seq[, i]
      prob <- 1 / n_row
      n_picked <- 0
      finished <- FALSE
      times_picked <- numeric(n_row)
      
      # Determine the number of times each row gets picked
      for (row in 1:n_row) {
        if (!finished) {
          possible_picks <- 0:(n_row - n_picked)
          pbinom_prob <- min(1, prob / (1 - (row - 1) * prob))
          y <- pbinom(possible_picks, size = n_row - n_picked,
                      prob = pbinom_prob)
          
          times_picked[row] <- which.max(y > qmc_seq[row]) - 1
          n_picked <- n_picked + times_picked[row]
          
          finished <- (n_picked == n_row)
        }
      }
      
      # Shuffle the picked rows
      ind <- unlist(mapply(rep, 1:n_row, times_picked)) 
      ind <- sample(ind, n_row)
    } else {
      ind <- sample(1:n_row, n_row, replace = TRUE)
    }
    
    # Map from original data to bootstrapped data
    orig_boots <- ind
    # Map from bootstrapped data back to original
    boots_orig <- match(1:length(ind), ind)
    
    # Store the mappings
    ind_boots[[paste0("bootstrap-", i)]] <- list(
      orig_boots = orig_boots,
      boots_orig = boots_orig
    )
  }
  
  output <- list(
    ind_boots = ind_boots,
    sobol_seq = sobol_seq
  )
  
  return(output)
}
