## Packages --------------------------------------------------------------------
# library(tidyverse)
# library(furrr)

## Functions to clean raw data -------------------------------------------------
# Normalize variables
normalize = function(data){
  
  out = data %>% 
    mutate_all(
      ~ (. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE)
    )
  
  out[is.na(out)] = 0
  return(out)
}

# Prepare auxiliary variables
prepare_auxiliary = function(auxiliary){
  
  if (!is.null(auxiliary)) {
    out = auxiliary %>% 
      normalize %>% 
      `colnames<-`(
        paste0(
          "auxiliary-",
          1:n_industry
        )
      )
  }else{
    out = NULL
  }
  
  return(out)
}

# Create polytime variables
create_polytime = function(n_polytime){
  
  polytime = 1:n_polytime %>% 
    as.list %>% 
    map(
      ~ (1:n_time)^.
    ) %>% 
    do.call("cbind", .) %>%
    data.frame %>% 
    mutate_all(
      ~ (. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE)
    ) %>% 
    `colnames<-`(
      paste0(
        "polytime-",
        1:n_polytime
      )
    )
  
  return(polytime)
}

# Clean data
clean_data = function(data){
  # scale down
  data$multiplier = data$detail %>% 
    colMeans(., na.rm = TRUE) %>% 
    mean(., na.rm = TRUE)
  
  data$detail = (data$detail/data$multiplier) %>% 
    `colnames<-`(
      paste0(
        "industry-",
        1:n_industry
      )
    )
  
  data$q_total = data$q_total/data$multiplier
  data$a_total = data$a_total/data$multiplier
  data$a_q_total = data$a_q_total/data$multiplier
  
  # add auxiliary variables
  if (!is.null(data$auxiliary)) {
    data$detail_with_aux = cbind(
      data$detail,
      prepare_auxiliary(data$auxiliary),
      create_polytime(n_polytime)
    )
  }else{
    data$detail_with_aux = cbind(
      data$detail,
      create_polytime(n_polytime)
    )
  }
  
  # add intercept
  data$detail_with_aux = data$detail_with_aux %>% 
    mutate(
      intercept = 1
    ) %>% 
    dplyr::select(intercept, everything())
  
  return(data)
}


## Functions to prepare Owen-scrambled sobol sequence --------------------------
prepare_sobol = function(dim, n = 1000, seed = 316){
  output = randtoolbox::sobol(
    n = n,
    dim = dim,
    init = TRUE,
    scrambling = 1, # Owen type
    seed = seed
  )
  return(output)
}


## Functions to do obs-by-obs bootstrapphing -----------------------------------
bootstrapping = function(nrows, n_bootstrap, sobol_seq = NULL, qmc = TRUE){
  # sobol seq
  if (is.null(sobol_seq)) {
    sobol_seq = prepare_sobol(
      dim = n_bootstrap,
      n = nrows,
      seed = as.integer(runif(1)*1000000)
    )
  }

  ind_boots = 1:n_bootstrap %>% 
    set_names(
      paste0("bootstrap-", 1:n_bootstrap)
    ) %>% 
    future_map(
      ~{
        if (qmc) {
          # start
          qmc_seq = sobol_seq[,.]
          prob = 1/nrows
          n_picked = 0
          finished = FALSE
          times_picked = rep(0, nrows)
          
          # determine # of times being picked
          for (row in 1:nrows) {
            if (!finished) {
              times_could_be_picked = 0:(nrows - n_picked)
              pbinom_prob = prob/(1-(row-1)*prob)
              if (pbinom_prob > 1) {
                pbinom_prob = 1
              }
              y = pbinom(
                q = times_could_be_picked,
                size = nrows - n_picked,
                prob = pbinom_prob
              )
              times_picked[row] = (y > qmc_seq[row]) %>% 
                which %>% 
                min - 1
              n_picked = n_picked + times_picked[row]
            }
            if (n_picked == nrows) {
              finished = TRUE
            }
          }
          
          # shuffle the picked rows
          ind = mapply(rep, 1:nrows, times_picked) %>% 
            do.call("c", .) %>% 
            sample(., nrows, replace = FALSE)
        }else{
          ind = sample(1:nrows, nrows, replace = TRUE)
        }
        
        # get indices for original to bootstrapped
        orig_boots = ind
        
        # get indices for boostrapped to original
        boots_orig = match(1:length(ind), ind)
        
        list(
          orig_boots = orig_boots,
          boots_orig = boots_orig
        )
      }
    )
  
  # finish
  output = list(
    ind_boots = ind_boots,
    sobol_seq = sobol_seq
  )
  
  return(output)
}


## Functions for E-step --------------------------------------------------------
# sweep sufficient statistics Q
sweep = function(Q, pos){
  nrows = nrow(Q)
  ncols = ncol(Q)
  if (nrows != ncols) {
    print("Error: input matrix must be symmetric.")
    return()
  }else{
    pos_compl = setdiff(1:nrows, pos)
    if (length(pos) == nrows) {
      result = try(
        {-solve(Q)},
        silent = TRUE
      )
      if ("try-error" %in% class(result)) {
        result = MASS::ginv(Q)
      }
      output = result
    }else if ((length(pos) == 1) & (1 %in% pos)) {
      Q = Q/Q[pos, pos]
      Q[pos_compl, pos_compl] = Q[pos_compl, pos_compl] - 
        Q[pos_compl, pos] %o% Q[pos, pos_compl]
      Q[pos, pos] = -1
      output = Q
    }else{
      result = try(
        {solve(Q[pos, pos])},
        silent = TRUE
      )
      if ("try-error" %in% class(result)) {
        result = MASS::ginv(Q[pos, pos])
      }
      Q[pos, pos] = result
      Q[pos, pos_compl] = Q[pos, pos] %*% Q[pos, pos_compl]
      Q[pos_compl, pos_compl] = Q[pos_compl, pos_compl] - 
        Q[pos_compl, pos] %*% Q[pos, pos_compl]
      Q[pos_compl, pos] = t(Q[pos, pos_compl])
      Q[pos, pos] = -Q[pos, pos]
      output = Q
    }
  }
  return(output)
}

# E-step
E_step = function(data){
  # add_var is used to store additional variance
  add_var = matrix(rep(0, nrow(data$Q)^2), nrow(data$Q), nrow(data$Q))
  
  for (type in 1:nrow(data$miss_types)) {
    # sweep Q for each miss type
    if (data$miss_types[type,] %>% as.integer %>% sum != 0) {
      ind = which(data$miss_type_ind == type)
      sweep_pos = which(data$miss_types[type,] == FALSE)[-1]
      theta = data$Q_old
      
      # sweep
      theta = sweep(theta, 1)
      theta = sweep(theta, sweep_pos) # try to feed auxiliary variables one by one
      
      # fill expectations
      observed = data$data_orig[ind,]*(!data$missing_orig[ind,])
      expectation = observed %*% theta
      data$data_orig[ind,][data$missing_orig[ind,]] = 
        expectation[data$missing_orig[ind,]]
      
      # store additional variance
      add_var = add_var + (t(matrix(data$missing_orig[ind,],nrow = length(ind))) %*%
                             matrix(data$missing_orig[ind,],nrow = length(ind))) * theta
    }
  }
  
  # construct Q_eStep
  Q_eStep = t(data$data_orig) %*% data$data_orig + add_var
  
  # get V matrix for multi-scale step
  data$V = diag(sweep(Q_eStep, 1))[2:(n_industry + 1)] %>% 
    rep(., each = 4) %>% 
    diag
  
  # finish
  return(data)
}


## Functions for multi-scale updating step -------------------------------------
# create transform matrix y_z
create_trans_y_z = function(n_time, n_year, n_industry){
  output = matrix(1:(n_time*n_industry), nrow = n_time) %>% 
    split(., rep(1:n_year, each = n_time/n_year)) %>% 
    lapply(., matrix, ncol = 1) %>% 
    do.call("cbind", .)
  return(output)
} 

# create transform matrix z_y
create_trans_z_y = function(n_time, n_year, n_industry){
  output = matrix(1:(n_time*n_industry), ncol = n_year) %>% 
    split(., rep(1:n_industry, each = n_time/n_year)) %>% 
    lapply(., matrix, ncol = 1) %>% 
    do.call("cbind", .)
  return(output)
}

# create transform matrix H
create_H = function(n_time, n_year, n_industry){
  n_time_per_year = n_time/n_year
  output = rbind(
    diag(1, n_time_per_year*n_industry),
    matlab::repmat(diag(1, n_time_per_year), 1, n_industry),
    kronecker(diag(1, n_industry), rep(1,4)) %>% t,
    rep(1, n_time_per_year*n_industry)
  )
  return(output)
}

# MSU_step
MSU_step = function(data){
  # transform data into z vector
  z = data$data_orig[, 2:(n_industry + 1)][data$y_z] %>% 
    matrix(., ncol = n_year) %>% 
    rbind(
      .,
      data$q_total,
      data$a_total,
      data$a_q_total
    )
  
  # store means and variances of the missing values
  data$gamma_m = NULL
  data$gamma_z = data$y_z * 0
  data$omega_m = NULL
  data$omega_y = NULL
  data$add_var_MSU = data$Q*0
  
  for (year in 1:n_year) {
    if (sum(data$missing_z[,year]) != 0) {
      # prepare z_target and H
      z_target = z[, year]
      H = data$H
      
      # delete missing totals
      available_z = (!is.na(z_target)) %>% which
      z_target = z_target[available_z]
      H = H[available_z,]
      
      # construct mu and sigma
      mu = H %*% z_target[1:(n_time_per_year*n_industry)]
      sigma = H %*% data$V %*% t(H)
      
      # get indices for missing and observed values
      ind_miss = c(
        data$missing_z[, year],
        rep(FALSE, length(z_target)-n_time_per_year*n_industry)
      )
      ind_obs = !ind_miss
      
      # partitions of sigma
      sigma_mm = sigma[ind_miss, ind_miss]
      sigma_oo = sigma[ind_obs, ind_obs]
      sigma_om = sigma[ind_obs, ind_miss]
      sigma_mo = t(sigma_om)
      
      # spectral decomposition
      # eig_result = matlib::Eigen(sigma_oo, max.iter = 100)    # may cause error
      eig_result = eigen(sigma_oo)
      sigma_oo_D = eig_result$values
      sigma_oo_P = eig_result$vectors
      
      
      non_zeros = sigma_oo_D > 0.0000001
      sigma_oo_D = sigma_oo_D[non_zeros]
      sigma_oo_P = sigma_oo_P[, non_zeros]
      
      # Moore-Penrose inverse
      sigma_oo_MPI = sigma_oo_P %*% diag(1/sigma_oo_D) %*% t(sigma_oo_P)
      
      # update the covariance matrix
      omega_m = sigma_mm - sigma_mo %*% sigma_oo_MPI %*% sigma_om
      
      # update the conditional mean
      gamma_m = mu[ind_miss] + 
        sigma_mo %*% sigma_oo_MPI %*% (z_target[ind_obs] - mu[ind_obs])
      
      # store gamma_m
      data$gamma_m[[year]] = gamma_m
      data$gamma_z[data$missing_z[, year], year] = gamma_m
      
      # store omega_m
      data$omega_m[[year]] = omega_m
      order_m = data$missing_orig[(n_time_per_year*(year-1)+1):(n_time_per_year*year), ]*1
      order_m[order_m == 1] = 1:sum(ind_miss)
      for (time in 1:n_time_per_year) {
        ind_Q = (order_m[time, ] != 0) %>% which
        ind_omega_m = order_m[time, ind_Q]
        data$omega_y[[n_time_per_year*(year-1)+time]] = data$Q*0
        data$omega_y[[n_time_per_year*(year-1)+time]][ind_Q, ind_Q] = 
          omega_m[ind_omega_m, ind_omega_m]
        data$add_var_MSU = data$add_var_MSU + data$omega_y[[n_time_per_year*(year-1)+time]]
      }
    }
  }
  
  # transform gamma to y format
  data$gamma_y = data$data_orig*0
  data$gamma_y[, 2:(n_industry+1)] = data$gamma_z[data$z_y %>% as.vector]
  
  return(data)
}


## Functions for M-step --------------------------------------------------------
M_step = function(data){
  # construct new sufficient statistics Q
  data$data_orig = data$data_orig*(!data$missing_orig) + data$gamma_y
  data$data_boots = data$data_orig[data$orig_boots,]
  data$Q = t(data$data_boots) %*% data$data_boots + data$add_var_MSU
  
  return(data)
}


## Functions to simulate missing values ----------------------------------------
rtmvn_sig = function(n, Mean, Sigma, D = diag(1, length(Mean)), lower = rep(0, length(Mean)), 
                     upper = rep(1000000, length(Mean)), int = Mean, burn = 10, thin = 1){
  if (length(Mean) == 1) {
    result <- tmvmixnorm::rtuvn(n = n, mean = Mean, sd = c(Sigma), lower = lower, 
                                upper = upper)
  }
  else {
    if (any(lower >= upper)) 
      stop("lower bound must be smaller than upper bound\n")
    # bound.check <- 0
    # if (!is.null(int)) {
    #   inits_test <- D %*% int
    #   lower.log <- inits_test >= lower + 1e-08
    #   upper.log <- inits_test <= upper - 1e-08
    #   bound.check <- prod(lower.log * upper.log)
    #   if (bound.check == 0) 
    #     cat("initial is outside or too close from boundary, will be auto-corrected by ginv()!\n")
    # }
    # else if (bound.check == 0) {
    #   D.inv <- MASS::ginv(D)
    #   int <- D.inv %*% (lower + upper)/2
    # }
    if (any(c(burn, thin, n)%%1 != 0)) 
      stop("burn, thin and n must be integer\n")
    if (any(c(burn, thin, n - 1) < 0)) 
      stop("burn, thin must be  non-negative interger, n must be positive integer\n")
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
          z[j] <- tmvmixnorm::rtuvn(lower = lower.j, upper = upper.j) # add QMC here
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

# simulate_missing_value_old = function(data, n_simulation, qmc = TRUE){
#   data$data_imputed = 1:n_simulation %>% 
#     as.list %>% 
#     set_names(paste0("simulation-", 1:n_simulation)) %>% 
#     map(
#       ~ data$data_orig*(!data$missing_orig)
#     )
#   for (year in 1:length(data$gamma_m)) {
#     if (length(data$gamma_m[[year]]) > 0) {
#       # prepare gamma_m and omega_m
#       gamma_m = data$gamma_m[[year]]
#       omega_m = data$omega_m[[year]]
#       
#       # spectral decomposition
#       eig_result = matlib::Eigen(omega_m, max.iter = 1000)    # may cause error
#       omega_m_D = eig_result$values
#       omega_m_P = eig_result$vectors
#       non_zeros = omega_m_D > 0.0000001
#       omega_m_D = omega_m_D[non_zeros]
#       omega_m_P = omega_m_P[, non_zeros]
#       
#       # random draws
#       if (qmc) {
#         random_seq = prepare_sobol(
#           dim = sum(non_zeros),
#           n = n_simulation,
#           seed = as.integer(runif(1)*1000000)
#         ) %>% 
#           matrix(., nrow = n_simulation)
#       }else{
#         random_seq = runif(sum(non_zeros) * n_simulation) %>% 
#           matrix(., nrow = n_simulation)
#       }
#       
#       # simulation
#       for (sim in 1:n_simulation) {
#         if (length(omega_m_D)==1) {
#           data$data_imputed[[sim]][data$y_z[, year][data$missing_z[, year]] + n_time] = 
#             gamma_m + omega_m_P * omega_m_D^0.5 * qnorm(random_seq[sim, ])
#         }else{
#           data$data_imputed[[sim]][data$y_z[, year][data$missing_z[, year]] + n_time] = 
#             gamma_m + omega_m_P %*% diag(omega_m_D)^0.5 %*% qnorm(random_seq[sim, ])
#         }
#       }
#     }
#   }
#   return(data)
# }

simulate_missing_value = function(data, n_simulation,
                                  lower = 0, upper = 100000,
                                  burn = 10, thin = 3){
  data$data_imputed = 1:n_simulation %>% 
    as.list %>% 
    set_names(paste0("simulation-", 1:n_simulation)) %>% 
    map(
      ~ data$data_orig*(!data$missing_orig)
    )
  for (year in 1:length(data$gamma_m)) {
    if (length(data$gamma_m[[year]]) > 0) {
      
      # prepare gamma_m and omega_m
      gamma_m = data$gamma_m[[year]]
      omega_m = data$omega_m[[year]]
      d = length(gamma_m)
      impute_ind = data$y_z[, year][data$missing_z[, year]] + n_time
      
      # simulation
      for (sim in 1:n_simulation) {
        ys = rtmvn_sig(n = 1,
                       Mean = gamma_m,
                       Sigma = omega_m,
                       lower = rep(lower, d),
                       upper = rep(upper, d),
                       burn = burn,
                       thin = thin)
        data$data_imputed[[sim]][impute_ind] = ys
      }
    }
  }
  return(data)
}


## MBEMMI ----------------------------------------------------------------------
MBEMMI = function(data_raw, n_time_per_year = 4, n_polytime = 3,
                  n_imputation = 10, n_simulation = 1,
                  qmc_boots = TRUE, qmc_simulation = TRUE,
                  tol = 0.00001, sobol_seq = NULL, max_inter = 10000){
  # parameters
  n_industry <<- data_raw$detail %>% ncol
  n_time <<- data_raw$detail %>% nrow
  n_time_per_year <<- n_time_per_year
  n_polytime <<- n_polytime
  n_imputation <<- n_imputation
  n_simulation <<- n_simulation
  n_year <<- n_time/n_time_per_year
  
  # clean data
  data_raw = clean_data(data_raw)
  
  # bootstrap indices
  boots_result = bootstrapping(n_time, n_imputation * 5,
                               sobol_seq = sobol_seq, qmc = qmc_boots) # get redundant ones for higher missing rate blocks
  boots_raw = boots_result$ind_boots
  boots_sobol_seq = boots_result$sobol_seq
  
  # missing types
  missing = is.na(data_raw$detail_with_aux)
  miss_types = unique(missing)
  miss_type_ind = prodlim::row.match(
    missing %>% data.frame,
    miss_types %>% data.frame
  )
  
  # filter out completely missing ind_boots
  ind_pass = NULL
  boots_pass = 0
  boots_count = 0
  while(boots_pass < n_imputation){
    boots_count = boots_count + 1
    check_missing = colSums(missing[boots_raw[[boots_count]]$orig_boots,])
    if (!any(check_missing == n_time)) {
      ind_pass = c(ind_pass, boots_count)
      boots_pass = boots_pass + 1
    }
  }
  ind_boots = boots_raw[ind_pass]
  
  # imputation
  result = 1:n_imputation %>% 
    as.list %>% 
    set_names(paste0("imputation-", 1:n_imputation)) %>% 
    future_map(
      ~{
        m = .
        ## start
        finished = FALSE
        converged = FALSE
        toasted = FALSE
        iter_count = 1
        
        # prepare data
        data = ind_boots[[m]]
        data$data_orig = data_raw$detail_with_aux %>% as.matrix
        data$data_boots = data_raw$detail_with_aux[data$orig_boots,] %>% as.matrix
        data$missing_orig = missing
        data$missing_boots = is.na(data$data_boots)
        data$data_orig[data$missing_orig] = 0
        data$data_boots[data$missing_boots] = 0
        data$miss_types = miss_types
        data$miss_type_ind = miss_type_ind
        data$y_z = create_trans_y_z(n_time, n_year, n_industry)
        data$z_y = create_trans_z_y(n_time, n_year, n_industry)
        data$missing_z = data$missing_orig[, 2:(n_industry + 1)][data$y_z] %>% 
          matrix(., ncol = n_year)
        data$q_total = data_raw$q_total[[1]] %>% matrix(., ncol = n_year)
        data$a_total = data_raw$a_total %>% t
        # data$a_q_total = data_raw$a_q_total[[1]] %>% t
        data$a_q_total = NULL
        data$H = create_H(n_time, n_year, n_industry)
        
        # fill the missing values with random draws
        data$data_boots[data$missing_boots] = mapply(
          rnorm,
          n_time,
          data$data_boots %>% 
            colMeans(., na.rm = TRUE),
          data$data_boots %>% 
            apply(., 2, sd, na.rm=TRUE)
        ) %>% 
          .[data$missing_boots]
        
        # construct initial sufficient statistics Q
        data$Q = t(data$data_boots) %*% data$data_boots
        data$Q_old = data$Q
        
        # start EMM
        cat(paste0("Imputation ", m, ": EMM algorithm start...\n"))
        while(!finished){
          # E-step
          data = E_step(data)
          
          # Multi-scale updating step
          data = MSU_step(data)
          
          # M-step
          data = M_step(data)
          
          # simulate missing values
          if (converged) {
            data = simulate_missing_value(data, n_simulation) # , qmc = qmc_simulation
            cat(paste0("Imputation ", m, ": EMM algorithm finished.\n"))
            finished = TRUE
          }else if (toasted){
            cat(paste0("Imputation ", m, ": EMM algorithm didn't converge.\n"))
            finished = TRUE
          }
          
          # Check convergence
          if (!finished) {
            Q_diff = abs(data$Q - data$Q_old) %>% 
              rowSums %>% 
              sum
            if (Q_diff < tol) {
              converged = TRUE
            }else if (Q_diff > 10000 | iter_count > max_inter){
              toasted = TRUE
            }
            # cat("-")
            print(paste0("Iteration: ", iter_count, ", Q_diff: ", Q_diff, "."))
            data$Q_old = data$Q
            iter_count = iter_count + 1
          }
        }
        
        # multipliers
        data$data_final = data$data_imputed %>% 
          map(
            ~ round(.[, 2:(n_industry + 1)]*data_raw$multiplier, 0)
          )
        
        data$converged = converged
        
        data
      }
    )
  
  return(result)
}

