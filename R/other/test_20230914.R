method = "MBEMMI"
n_imp = 1
parallel_imp = FALSE
parallel_level = FALSE

data_qcew = list_data_qcew[[8]]

### preparation

# define map type
fun_map_imp <- ifelse(parallel_imp, future_map, map)
fun_map_level <- ifelse(parallel_level, future_map2, map2)

# partition data
blocks <- partition_data(data_qcew)

# generate sobol sequence
n_time <- nrow(blocks[[1]][[1]][[1]])
n_block <- sum(sapply(blocks, length))
sobol_seq <- prepare_sobol(n = n_time * n_block, dim = n_imp * 5)
list_sobol_imp <- lapply(1:n_imp, function(i){
  start_col <- (i - 1) * 5 + 1
  end_col <- min(i * 5, n_imp * 5)
  sobol_seq[, start_col:end_col, drop = FALSE]
})







# parallel imputations


sobol_imp <- list_sobol_imp[[1]]
data_imputed <- data_qcew$suppressed$qtrly

# sequential levels
for (level in 1:length(blocks)) {
  list_data_block <- blocks[[level]]
  n_block_level <- length(list_data_block)
  if (n_block_level == 0) next
  
  # prepare sobol sequences for blocks
  list_sobol_block <- lapply(1:n_block_level, function(i){
    start_row <- (i - 1) * n_time + 1
    end_row <- min(i * n_time, n_block_level * n_time)
    sobol_imp[start_row:end_row, , drop = FALSE]
  })
  sobol_imp <- sobol_imp[-(1:n_block_level * n_time), ]
  
  # parallel blocks
  result_level <- NULL
  
  # Loop over each name in data_level
  for (item in names(list_data_block)) {
    ind = which(names(list_data_block) == item)
    data <- list_data_block[[item]]
    sobol_block <- list_sobol_block[[ind]]
    code_higher <- colnames(data$q_total)
    code_lower <- colnames(data$detail)
    
    # fill missing q_total
    if (any(is.na(data$q_total))) {  
      data$q_total <- data_imputed %>%
        filter(industry_code %in% code_higher) %>%
        dplyr::select(-c("industry_code", "industry_title")) %>%
        t() %>%
        as.data.frame() %>%
        dplyr::rename_all(~code_higher)
    }
    
    # impute missing values
    if (length(code_lower) > 1) {
      if (method == "MBEMMI") {
        res <- MBEMMI(data_raw = data, 
                      n_imp = 1,
                      sobol_seq = sobol_block)[[1]]
        # if didn't converge, use new sobol_seq
        while (res$converged == FALSE) {
          res <- MBEMMI(data_raw = data, 
                        n_imp = 1,
                        sobol_seq = NULL)[[1]]
        }
        result_level[[item]] <- res[["data_final"]][[1]]
      }else if(method == "BMMI"){
        result_level[[item]] <- BMMI(data, n_imputation = 1)[[1]]
      }else if(method == "EMB"){
        res <- EMB(data, n_imputation = 1)[[1]]
        if (length(res) > 0) {
          result_level[[item]] <- res
        }else{
          result_level[[item]] <- data$detail
        }
      }
    }else{
      result_level[[item]] <- list()
    }
    
    # Extract column names for 'detail' and 'q_total'
    code_lower <- colnames(list_data_block[[item]]$detail)
    code_higher <- colnames(list_data_block[[item]]$q_total)
    
    # Identify rows in data_imputed that match code_lower/code_higher
    ind_lower <- which(data_imputed$industry_code %in% code_lower)
    ind_higher <- which(data_imputed$industry_code %in% code_higher)
    
    # Depending on the condition, update these rows
    if (length(code_lower) > 1) {
      data_imputed[ind_lower, 3:ncol(data_imputed)] <-
        t(result_level[[item]])
    } else {
      data_higher <- data_imputed[ind_higher, 3:ncol(data_imputed)]
      data_imputed[ind_lower, 3:ncol(data_imputed)] <- data_higher
    }
  }
}

for (i in 1:1000) {
  print(i)
  set.seed(i)
  res <- MBEMMI(data_raw = data, 
                n_imp = 1,
                sobol_seq = NULL)[[1]]
}

# level = 4
# item = "22111"
# i = 7


# data_error = data
set.seed(i)
data_raw = data_error
n_time_per_year = 4
n_polytime = 3
n_imp = 1
n_sim = 1
qmc_boots = TRUE
sobol_seq = sobol_block
sobol_max_burn = 1000
sobol_scrambling = 0
sobol_seed = 316
parallel = FALSE
qmc_simulation = TRUE
tol = 0.00001
max_inter = 10000




# Define parameters
n_industry <- ncol(data_raw$detail)
n_time <- nrow(data_raw$detail)

# Clean the raw data
data_raw <- clean_data(data_raw, n_polytime = n_polytime)

# Generate bootstrapping indices
boots_result <- bootstrapping(
  n_row = n_time,
  n_bootstrap = n_imp * 5,
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
boots_ind <- boots_ind[!miss_all][1:n_imp]


# Imputation
m <- 1
boots <- boots_ind[[1]]
finished <- FALSE
converged <- FALSE
toasted <- FALSE
iter_count <- 1

# Prepare data
data <- prepare_data(data_raw, boots, miss_bool,
                     miss_types, miss_type_ind,
                     n_time, n_time_per_year, n_industry)

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
      toasted <- TRUE
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





a = sim_miss(data, n_sim) 



# data0 = data
data = data0
lower = 0
upper = 100000
burn = 10
thin = 3

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
  impute_ind <- data$y_z[, year][data$miss_z[, year]] + data$n_time
  
  # simulation
  for (sim in 1:n_sim) {
    ys <- rtmvn_sig(
      n = 1,
      Mean = gamma_m,
      Sigma = round(omega_m, 10),
      lower = rep(lower, d),
      upper = rep(upper, d),
      burn = burn,
      thin = thin
    )
    data$data_imputed[[sim]][impute_ind] <- ys
  }
}

ys <- rtmvn_sig(
  n = 1,
  Mean = gamma_m,
  Sigma = round(omega_m, 10),
  lower = rep(lower, d),
  upper = rep(upper, d),
  burn = burn,
  thin = thin
)

n = 1
Mean = gamma_m
Sigma = round(omega_m, 10)
lower = rep(lower, d)
upper = rep(upper, d)
burn = burn
thin = thin
D = diag(1, length(Mean))
int = Mean

