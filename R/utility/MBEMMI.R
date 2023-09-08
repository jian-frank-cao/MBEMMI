## Functions to clean raw data -------------------------------------------------
# Function to normalize the data
normalize <- function(data) {
  
  # Subtract mean and divide by standard deviation for each column
  normalized_data <- data %>%
    mutate(across(everything(), 
                  ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)
    ))
  
  # Replace NA values with 0                    ### NA may be caused by sd = 0  
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
prepare_sobol <- function(dim, n = 1000, max_burn = 1000, 
                          scrambling = 0, seed = 316){
  
  # Randomly decide on a 'burn-in' period between 0 and max_burn.
  burn_in <- sample(0:max_burn, 1)
  
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


## Functions to perform observation-by-observation bootstrapping ---------------
bootstrapping <- function(n_row, n_bootstrap, qmc = TRUE,
                          sobol_seq = NULL, max_burn = 1000,
                          scrambling = 0, seed = 316){
  
  # If no Sobol sequence is provided, generate one
  if (is.null(sobol_seq)) {
    sobol_seq <- prepare_sobol(
      dim = n_bootstrap,
      n = n_row,
      max_burn = max_burn,
      scrambling = scrambling,
      seed = seed
    )
  }
  
  boots_ind <- list()
  
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
    boots_orig <- match(1:length(ind), ind)             ### remember duplicates
    
    # Store the mappings
    boots_ind[[paste0("bootstrap-", i)]] <- list(
      orig_boots = orig_boots,
      boots_orig = boots_orig
    )
  }
  
  output <- list(
    boots_ind = boots_ind,
    sobol_seq = sobol_seq
  )
  
  return(output)
}


## Functions to prepare data for imputation ------------------------------------
prepare_data <- function(data_raw, boots, miss_bool, miss_types,
                         miss_type_ind, n_time, n_year, n_industry) {
  
  # Initial data assignments
  data_orig <- as.matrix(data_raw$detail_with_aux)
  data_boots <- as.matrix(data_raw$detail_with_aux[boots$orig_boots,])
  miss_orig <- miss_bool
  miss_boots <- is.na(data_boots)
  
  # Setting missing values to zero
  data_orig[miss_orig] <- 0
  data_boots[miss_boots] <- 0
  
  # Compute means and standard deviations for bootstrapped data
  data_boots_means <- colMeans(data_boots, na.rm = TRUE)
  data_boots_sds <- apply(data_boots, 2, sd, na.rm=TRUE)
  
  # Fill the missing values in bootstrapped data with random draws
  data_boots[miss_boots] <- mapply(rnorm, n_time, data_boots_means,
                                   data_boots_sds)[miss_boots]
  
  # Obtain transformation indices y_z and z_y
  y_z <- create_trans_y_z(n_time, n_year, n_industry)
  z_y <- create_trans_z_y(n_time, n_year, n_industry)
  
  # Construct the main data list
  data <- list(
    orig_boots = boots$orig_boots,
    boots_orig = boots$boots_orig,
    data_orig = data_orig,
    data_boots = data_boots,
    miss_orig = miss_orig,
    miss_boots = miss_boots,
    miss_types = miss_types,
    miss_type_ind = miss_type_ind,
    y_z = y_z,
    z_y = z_y,
    miss_z = matrix(miss_orig[, 2:(n_industry + 1)][y_z], ncol = n_year),
    q_total = matrix(data_raw$q_total[[1]], ncol = n_year),
    a_total = t(data_raw$a_total),
    a_q_total = NULL,
    H = create_H(n_time, n_year, n_industry)
  )
  
  return(data)
}


## Functions for E-step --------------------------------------------------------
# sweep sufficient statistics Q
sweep <- function(Q, pos) {
  # Ensure the matrix is square
  if (nrow(Q) != ncol(Q)) {
    stop("Error: input matrix must be symmetric.")
  }
  
  nrows <- nrow(Q)
  pos_compl <- setdiff(1:nrows, pos)
  
  # Case 1: Length of position is equal to number of rows
  if (length(pos) == nrows) {
    result <- tryCatch(
      -solve(Q),
      error = function(e) MASS::ginv(Q)
    )
    return(result)
  }
  
  # Case 2: Single position, and it's the first one
  if ((length(pos) == 1) && (1 %in% pos)) {
    Q <- Q / Q[pos, pos]
    Q[pos_compl, pos_compl] <- Q[pos_compl, pos_compl] -
      Q[pos_compl, pos] %o% Q[pos, pos_compl]
    Q[pos, pos] <- -1
    return(Q)
  }
  
  # Case 3: General scenario
  result <- tryCatch(
    solve(Q[pos, pos]),
    error = function(e) MASS::ginv(Q[pos, pos])
  )
  
  Q[pos, pos] <- result
  Q[pos, pos_compl] <- Q[pos, pos] %*% Q[pos, pos_compl]
  Q[pos_compl, pos_compl] <- Q[pos_compl, pos_compl] -
    Q[pos_compl, pos] %*% Q[pos, pos_compl]
  Q[pos_compl, pos] <- t(Q[pos, pos_compl])
  Q[pos, pos] <- -Q[pos, pos]
  
  return(Q)
}

# E-step
E_step <- function(data) {
  
  # Initialize matrix to store additional variance
  nQ <- nrow(data$Q)
  add_var <- matrix(0, nQ, nQ)
  
  # Loop through each missing type
  for (type in 1:nrow(data$miss_types)) {
    
    # If current missing type has any missing values
    if (sum(as.integer(data$miss_types[type, ])) != 0) {
      ind <- which(data$miss_type_ind == type)
      sweep_pos <- which(!data$miss_types[type, ])[-1]
      theta <- data$Q_old
      
      # Sweep operations
      theta <- sweep(theta, 1)
      theta <- sweep(theta, sweep_pos)   ### try auxiliary variables one by one
      
      # Compute expectations
      observed <- data$data_orig[ind, ] * !data$miss_orig[ind, ]
      expectation <- observed %*% theta
      data$data_orig[ind, ][data$miss_orig[ind, ]] <-
        expectation[data$miss_orig[ind, ]]
      
      # Update additional variance
      miss_mat <- matrix(data$miss_orig[ind, ], nrow = length(ind))
      add_var <- add_var + (t(miss_mat) %*% miss_mat) * theta
    }
  }
  
  # Compute Q matrix for E step
  Q_eStep <- t(data$data_orig) %*% data$data_orig + add_var
  
  # Compute V matrix for multi-scale step
  data$V <- diag(sweep(Q_eStep, 1))[2:(n_industry + 1)] %>% 
    rep(., each = 4) %>% 
    diag
  
  return(data)
}


## Functions for multi-scale updating step -------------------------------------
# create transform matrix y_z
create_trans_y_z <- function(n_time, n_year, n_industry){
  
  output <- matrix(1:(n_time * n_industry), nrow = n_time) %>% 
    split(., rep(1:n_year, each = n_time/n_year)) %>% 
    lapply(., matrix, ncol = 1) %>% 
    do.call("cbind", .)
  
  return(output)
} 

# create transform matrix z_y
create_trans_z_y <- function(n_time, n_year, n_industry){
  
  output <- matrix(1:(n_time * n_industry), ncol = n_year) %>% 
    split(., rep(1:n_industry, each = n_time/n_year)) %>% 
    lapply(., matrix, ncol = 1) %>% 
    do.call("cbind", .)
  
  return(output)
}

# create transform matrix H
create_H <- function(n_time, n_year, n_industry){
  n_time_per_year <- n_time/n_year
  output <- rbind(
    diag(1, n_time_per_year * n_industry),
    matlab::repmat(diag(1, n_time_per_year), 1, n_industry),
    kronecker(diag(1, n_industry), rep(1,4)) %>% t,
    rep(1, n_time_per_year * n_industry)
  )
  
  return(output)
}

# MSU_step
MSU_step <- function(data){
  # transform data into z vector
  z <- data$data_orig[, 2:(n_industry + 1)][data$y_z] %>% 
    matrix(., ncol = n_year) %>% 
    rbind(., data$q_total, data$a_total, data$a_q_total)
  
  # store means and variances of the missing values
  data$gamma_m <- NULL
  data$gamma_z <- data$y_z * 0
  data$omega_m <- NULL
  data$omega_y <- NULL
  data$add_var_MSU <- data$Q * 0
  
  for (year in 1:n_year) {
    if (sum(data$miss_z[, year]) != 0) {
      # prepare z_target and H
      z_target <- z[, year]
      H <- data$H
      
      # delete missing totals
      available_z <- (!is.na(z_target)) %>% which
      z_target <- z_target[available_z]
      H <- H[available_z,]
      
      # construct mu and sigma
      mu <- H %*% z_target[1:(n_time_per_year*n_industry)]
      sigma <- H %*% data$V %*% t(H)
      
      # get indices for missing and observed values
      ind_miss <- c(
        data$miss_z[, year],
        rep(FALSE, length(z_target) - n_time_per_year * n_industry)
      )
      ind_obs <- !ind_miss
      
      # partitions of sigma
      sigma_mm <- sigma[ind_miss, ind_miss]
      sigma_oo <- sigma[ind_obs, ind_obs]
      sigma_om <- sigma[ind_obs, ind_miss]
      sigma_mo <- t(sigma_om)
      
      # spectral decomposition
      # eig_result = matlib::Eigen(sigma_oo, max.iter = 100)  # may cause error
      eig_result <- eigen(sigma_oo)
      sigma_oo_D <- eig_result$values
      sigma_oo_P <- eig_result$vectors
      
      non_zeros <- sigma_oo_D > 0.0000001
      sigma_oo_D <- sigma_oo_D[non_zeros]
      sigma_oo_P <- sigma_oo_P[, non_zeros]
      
      # Moore-Penrose inverse
      sigma_oo_MPI <- sigma_oo_P %*% diag(1/sigma_oo_D) %*% t(sigma_oo_P)
      
      # update the covariance matrix
      omega_m <- sigma_mm - sigma_mo %*% sigma_oo_MPI %*% sigma_om
      
      # update the conditional mean
      gamma_m <- mu[ind_miss] + 
        sigma_mo %*% sigma_oo_MPI %*% (z_target[ind_obs] - mu[ind_obs])
      
      # store gamma_m
      data$gamma_m[[year]] <- gamma_m
      data$gamma_z[data$miss_z[, year], year] <- gamma_m
      
      # store omega_m
      data$omega_m[[year]] <- omega_m
      order_m_start <- (n_time_per_year * (year - 1) + 1)
      order_m_end <- (n_time_per_year * year)
      order_m <- data$miss_orig[order_m_start:order_m_end, ] * 1
      order_m[order_m == 1] <- 1:sum(ind_miss)
      
      for (time in 1:n_time_per_year) {
        ind_Q <- (order_m[time, ] != 0) %>% which
        ind_omega_m <- order_m[time, ind_Q]
        
        data$omega_y[[n_time_per_year * (year - 1) + time]] <- data$Q * 0
        data$omega_y[[n_time_per_year * (year - 1) + time]][ind_Q, ind_Q] <- 
          omega_m[ind_omega_m, ind_omega_m]
        data$add_var_MSU <- data$add_var_MSU + 
          data$omega_y[[n_time_per_year * (year - 1) + time]]
      }
    }
  }
  
  # transform gamma to y format
  data$gamma_y <- data$data_orig * 0
  data$gamma_y[, 2:(n_industry + 1)] <- data$gamma_z[as.vector(data$z_y)]
  
  return(data)
}


## Functions for M-step --------------------------------------------------------
M_step <- function(data){
  
  # construct new sufficient statistics Q
  data$data_orig <- data$data_orig * (!data$miss_orig) + data$gamma_y
  data$data_boots <- data$data_orig[data$orig_boots, ]
  data$Q <- t(data$data_boots) %*% data$data_boots + data$add_var_MSU
  
  return(data)
}


## Functions to simulate missing values ----------------------------------------
# Modify rtmvn_sig function
rtmvn_sig <- function(n, Mean, Sigma, D = diag(1, length(Mean)), 
                     lower = rep(0, length(Mean)), 
                     upper = rep(1000000, length(Mean)), 
                     int = Mean, burn = 10, thin = 1){
  if (length(Mean) == 1) {
    result <- tmvmixnorm::rtuvn(n = n, mean = Mean, sd = c(Sigma),
                                lower = lower, upper = upper)
  }
  else {
    if (any(lower >= upper)) 
      stop("lower bound must be smaller than upper bound\n")
    if (any(c(burn, thin, n)%%1 != 0)) 
      stop("burn, thin and n must be integer\n")
    if (any(c(burn, thin, n - 1) < 0)) 
      stop("burn, thin must be  non-negative interger,
           n must be positive integer\n")
    if (is.vector(D) == TRUE) {
      Rtilde <- t(as.matrix(D))
      lower <- as.vector(lower)
      upper <- as.vector(upper)
    }else{
      Rtilde <- D
    }
    a <- lower - Rtilde %*% Mean
    b <- upper - Rtilde %*% Mean
    # gchol -------->
    eig_result = eigen(Sigma)
    Sigma_D = round(eig_result$values, 12)
    Sigma_P = round(eig_result$vectors, 12)
    Sigma.chol = Sigma_P %*% diag(Sigma_D)^0.5
    # end ----------<
    R <- Rtilde %*% Sigma.chol
    p <- ncol(R)
    # use zeros as initials ---------->
    z <- rep(0, length(Mean))
    # end ----------------------------<
    keep.x <- matrix(NA, ncol = p, nrow = (thin + 1) * n + burn)
    for (i in 1:((thin + 1) * n + burn)) {
      for (j in 1:p) {
        rj <- as.vector(R[, j])
        Rj <- as.matrix(R[, -j])
        zj <- as.vector(z[-j])
        a.temp <- a - Rj %*% zj
        b.temp <- b - Rj %*% zj
        pos <- rj > 0
        neg <- rj < 0
        if (sum(pos) == 0) {
          lower.pos <- -Inf
          upper.pos <- Inf
        }
        else {
          lower.pos <- max(a.temp[pos]/rj[pos])
          upper.pos <- min(b.temp[pos]/rj[pos])
        }
        if (sum(neg) == 0) {
          upper.neg <- Inf
          lower.neg <- -Inf
        }
        else {
          upper.neg <- min(a.temp[neg]/rj[neg])
          lower.neg <- max(b.temp[neg]/rj[neg])
        }
        lower.j <- max(lower.pos, lower.neg)
        upper.j <- min(upper.pos, upper.neg)
        if (lower.j <= upper.j) {
          z[j] <- tmvmixnorm::rtuvn(lower = lower.j, upper = upper.j) # add QMC
        }else{
          z[j] = lower.j
        }
      }
      x <- Sigma.chol %*% z + Mean
      keep.x[i, ] <- x
    }
    final.ind <- 1:((thin + 1) * n + burn)
    final.ind <- final.ind[(burn + 1):length(final.ind)]
    final.ind <- seq(1, length(final.ind), by = thin + 1) + 
      thin + burn
    if (n == 1) {
      result <- c(keep.x[final.ind, ])
    }
    else {
      result <- keep.x[final.ind, ]
    }
  }
  return(result)
}

# Simulate missing values                                             # add QMC
sim_miss <- function(data, n_sim,
                    lower = 0, upper = 100000,
                    burn = 10, thin = 3){
  data$data_imputed <- lapply(
    1:n_sim, 
    function(x) data$data_orig * (!data$miss_orig)
  )
  names(data$data_imputed) <- paste0("simulation-", 1:n_sim)
  
  for (year in 1:length(data$gamma_m)) {
    if (length(data$gamma_m[[year]]) == 0) next
    
    # prepare gamma_m and omega_m
    gamma_m <- data$gamma_m[[year]]
    omega_m <- data$omega_m[[year]]
    d <- length(gamma_m)
    impute_ind <- data$y_z[, year][data$miss_z[, year]] + n_time
    
    # simulation
    for (sim in 1:n_sim) {
      ys <- rtmvn_sig(
        n = 1,
        Mean = gamma_m,
        Sigma = omega_m,
        lower = rep(lower, d),
        upper = rep(upper, d),
        burn = burn,
        thin = thin
      )
      data$data_imputed[[sim]][impute_ind] <- ys
    }
  }
  
  return(data)
}


## MBEMMI ----------------------------------------------------------------------
MBEMMI <- function(data_raw, n_time_per_year = 4, n_polytime = 3,
                   n_imputation = 10, n_sim = 1,
                   qmc_boots = TRUE, sobol_seq = NULL, 
                   sobol_max_burn = 1000, sobol_scrambling = 0,
                   sobol_seed = 316, parallel = FALSE,
                   qmc_simulation = TRUE, tol = 0.00001, max_inter = 10000){
  
  # Define parameters
  n_industry <- ncol(data_raw$detail)
  n_time <- nrow(data_raw$detail)
  
  # Clean the raw data
  data_raw <- clean_data(data_raw)
  
  # Generate bootstrapping indices
  boots_result <- bootstrapping(
    n_row = n_time,
    n_bootstrap = n_imputation * 5,
    qmc = qmc_boots,
    sobol_seq = sobol_seq,
    max_burn = sobol_max_burn,
    scrambling = sobol_scrambling,
    seed = sobol_seed
  )
  boots_ind <- boots_result$boots_ind

  # Determine missing data types
  miss_bool <- is.na(data_raw$detail_with_aux)
  miss_types <- unique(miss_bool)
  miss_type_ind <- prodlim::row.match(as.data.frame(miss_bool),
                                      as.data.frame(miss_types))
  
  # Filter out bootstraps with complete missing data
  miss_all <- sapply(boots_ind, function(boot) {
    any(colSums(miss_bool[boot$orig_boots, ]) == n_time)
  })
  boots_ind <- boots_ind[!miss_all][1:n_imputation]
  
  
  # Imputation
  fun_map <- ifelse(parallel, future_map2, map2)
  result <- fun_map(
    as.list(1:n_imputation),
    boots_ind,
    ~{
      # Define parameters
      m <- .x
      boots <- .y
      finished <- FALSE
      converged <- FALSE
      toasted <- FALSE
      iter_count <- 1
      
      # Prepare data
      data <- prepare_data(data_raw, boots, miss_bool,
                           miss_types, miss_type_ind,
                           n_time, n_year, n_industry)
      
      # Construct initial sufficient statistics Q
      data$Q <- t(data$data_boots) %*% data$data_boots
      data$Q_old <- data$Q
      
      # Start EMM algorithm
      cat(paste0("Imputation ", m, ": EMM algorithm start...\n"))
      while(!finished){
        # E-step
        data <- E_step(data)
        
        # Multi-scale updating step
        data <- MSU_step(data)
        
        # M-step
        data <- M_step(data)
        
        # simulate missing values
        if (converged) {
          data <- sim_miss(data, n_sim)       ### qmc = qmc_simulation
          cat(paste0("Imputation ", m, ": EMM algorithm finished.\n"))
          finished <- TRUE
        }else if (toasted){
          cat(paste0("Imputation ", m, ": EMM algorithm didn't converge.\n"))
          finished <- TRUE
        }
        
        # Check convergence
        if (!finished) {
          Q_diff <- sum(rowSums(abs(data$Q - data$Q_old)))
          
          if (Q_diff < tol) {
            converged <- TRUE
          }else if (Q_diff > 10000 | iter_count > max_inter){
            toasted = TRUE
          }

          print(paste0("Iteration: ", iter_count, ", Q_diff: ", Q_diff, "."))
          data$Q_old <- data$Q
          iter_count <- iter_count + 1
        }
      }
      
      # Rescale data
      data$data_final <- lapply(data$data_imputed, function(x) {
        round(x[, 2:(n_industry + 1)] * data_raw$multiplier, 0)
      })
      
      # Finish
      data$converged <- converged
      data
    }
  )
  
  return(result)
}





