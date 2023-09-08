# partition data
blocks = partition_data(data_qcew)


## preparation -----------------------------------------------------------------

method = "MBEMMI"
imputation = NULL
n_imputation = 1
level = 2
item = "212"
data_level = blocks[[level]]
data = data_level[[item]]
code_higher = colnames(data$q_total)
code_lower = colnames(data$detail)
# result = MBEMMI(data, n_imputation = n_imputation)[[1]][["data_final"]][[1]]

data_raw = data
n_time_per_year = 4
n_polytime = 3
n_imputation = 10
n_sim = 1
qmc_boots = TRUE
sobol_seq = NULL
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


# Define parameters
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


## Start -----------------------------------------------------------------------

# Start EMM algorithm
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




## sim_miss --------------------------------------------------------------------
# data_error = data
data = data_error
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
      Sigma = omega_m,
      lower = rep(lower, d),
      upper = rep(upper, d),
      burn = burn,
      thin = thin
    )
    data$data_imputed[[sim]][impute_ind] <- ys
  }
}


## rtmvn_sig -------------------------------------------------------------------
n = 1
Mean = gamma_m
Sigma = round(omega_m, 15)
D = diag(1, length(Mean))
lower = rep(lower, d)
upper = rep(upper, d)
int = Mean
burn = burn 
thin = thin


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

