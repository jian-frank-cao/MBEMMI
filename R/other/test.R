# partition data
blocks = partition_data(data_qcew)

# imputation
result = 1:n_imputation %>% 
  as.list %>% 
  set_names(paste0("imputation-", 1:n_imputation)) %>% 
  future_map(
    ~{
      data_imputed = data_qcew$suppressed$qtrly
      for (level in 1:length(blocks)) {
        data_level = blocks[[level]]
        if (length(data_level) > 0) {
          result = data_level %>% 
            future_map(
              ~{
                data = .
                code_higher = colnames(data$q_total)
                code_lower = colnames(data$detail)
                if (any(is.na(data$q_total))) {  
                  data$q_total = data_imputed %>%
                    filter(industry_code %in% code_higher) %>% 
                    dplyr::select(-c("industry_code", "industry_title")) %>% 
                    t %>% 
                    data.frame %>% 
                    `colnames<-`(code_higher)
                }
                if (length(code_lower) > 1) {
                  if (method == "MBEMMI") {
                    MBEMMI(data, n_imputation = 1)[[1]][["data_final"]][[1]] # solve sobol seq problem + add qmc at simulation
                  }else if(method == "BMMI"){
                    BMMI(data, n_imputation = 1)[[1]]
                  }else if(method == "EMB"){
                    EMB(data, n_imputation = 1)[[1]]
                  }
                }else{
                  list()
                }
              }
            )
          for (item in names(data_level)) {
            code_lower = colnames(data_level[[item]]$detail)
            code_higher = colnames(data_level[[item]]$q_total)
            if (length(code_lower) > 1) {
              data_imputed[data_imputed$industry_code %in% code_lower,
                           3:ncol(data_imputed)] = result[[item]] %>% t
            }else{
              data_imputed[data_imputed$industry_code %in% code_lower,
                           3:ncol(data_imputed)] = data_imputed[data_imputed$industry_code %in% code_higher,
                                                                3:ncol(data_imputed)]
            }
          }
        }
      }
      data_imputed
    }
  )



method = "MBEMMI"
imputation = NULL
n_imputation = 1

data_imputed = data_qcew$suppressed$qtrly
for (level in 1:length(blocks)) {
  data_level = blocks[[level]]
  if (length(data_level) > 0) {
    for (item in names(data_level)) {
      data = data_level[[item]]
      code_higher = colnames(data$q_total)
      code_lower = colnames(data$detail)
      if (any(is.na(data$q_total))) {  
        data$q_total = data_imputed %>%
          filter(industry_code %in% code_higher) %>% 
          dplyr::select(-c("industry_code", "industry_title")) %>% 
          t %>% 
          data.frame %>% 
          `colnames<-`(code_higher)
      }
      if (length(code_lower) > 1) {
        if (method == "MBEMMI") {
          result = MBEMMI(data, n_imputation = 1)[[1]][["data_final"]][[1]] # solve sobol seq problem + add qmc at simulation
        }else if(method == "BMMI"){
          result = BMMI(data, n_imputation = 1)[[1]]
        }else if(method == "EMB"){
          result = EMB(data, n_imputation = 1)[[1]]
        }
      }else{
        result = list()
      }
    }
    
    # result = data_level %>% 
    #   future_map(
    #     ~{
    #       data = .
    #       code_higher = colnames(data$q_total)
    #       code_lower = colnames(data$detail)
    #       if (any(is.na(data$q_total))) {  
    #         data$q_total = data_imputed %>%
    #           filter(industry_code %in% code_higher) %>% 
    #           dplyr::select(-c("industry_code", "industry_title")) %>% 
    #           t %>% 
    #           data.frame %>% 
    #           `colnames<-`(code_higher)
    #       }
    #       if (length(code_lower) > 1) {
    #         if (method == "MBEMMI") {
    #           MBEMMI(data, n_imputation = 1)[[1]][["data_final"]][[1]] # solve sobol seq problem + add qmc at simulation
    #         }else if(method == "BMMI"){
    #           BMMI(data, n_imputation = 1)[[1]]
    #         }else if(method == "EMB"){
    #           EMB(data, n_imputation = 1)[[1]]
    #         }
    #       }else{
    #         list()
    #       }
    #     }
    #   )
    
    
    # for (item in names(data_level)) {
    #   code_lower = colnames(data_level[[item]]$detail)
    #   code_higher = colnames(data_level[[item]]$q_total)
    #   if (length(code_lower) > 1) {
    #     data_imputed[data_imputed$industry_code %in% code_lower,
    #                  3:ncol(data_imputed)] = result[[item]] %>% t
    #   }else{
    #     data_imputed[data_imputed$industry_code %in% code_lower,
    #                  3:ncol(data_imputed)] = data_imputed[data_imputed$industry_code %in% code_higher,
    #                                                       3:ncol(data_imputed)]
    #   }
    # }
  }
}


level = 2
item = "212"
data_level = blocks[[level]]
data = data_level[[item]]
code_higher = colnames(data$q_total)
code_lower = colnames(data$detail)
if (any(is.na(data$q_total))) {  
  data$q_total = data_imputed %>%
    filter(industry_code %in% code_higher) %>% 
    dplyr::select(-c("industry_code", "industry_title")) %>% 
    t %>% 
    data.frame %>% 
    `colnames<-`(code_higher)
}
if (length(code_lower) > 1) {
  if (method == "MBEMMI") {
    result = MBEMMI(data, n_imputation = 1)[[1]][["data_final"]][[1]] # solve sobol seq problem + add qmc at simulation
  }else if(method == "BMMI"){
    result = BMMI(data, n_imputation = 1)[[1]]
  }else if(method == "EMB"){
    result = EMB(data, n_imputation = 1)[[1]]
  }
}else{
  result = list()
}




output = randtoolbox::sobol(
  n = 100,
  dim = 3,
  init = TRUE,
  scrambling = 2, # Owen type
  seed = 123
)




level = 2
item = "212"
data_level = blocks[[level]]
data = data_level[[item]]
data2= data
data_raw <- data

n_industry <<- data_raw$detail %>% ncol
n_time <<- data_raw$detail %>% nrow
n_time_per_year <<- 4
n_polytime <<- 3
n_imputation <<- 10
n_simulation <<- 1
n_year <<- n_time/n_time_per_year

data_raw = clean_data(data_raw)

qmc_boots = TRUE
qmc_simulation = TRUE
tol = 0.00001
sobol_seq = NULL          ### change this
max_inter = 10000
sobol_max_burn = 1000
sobol_scrambling = 0
sobol_seed = 316

n_row = n_time,
n_bootstrap = n_imputation * 5,
qmc = qmc_boots,
sobol_seq = sobol_seq,
max_burn = sobol_max_burn,
scrambling = sobol_scrambling,
seed = sobol_seed

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

a = sweep(theta, 1)
b = sweep_old(theta, 1)
sum(colSums(a == b))
theta = sweep(theta, sweep_pos) 
a = sweep(theta, sweep_pos)
b = sweep_old(theta, sweep_pos)
sum(colSums(a == b))


data = data0
a = E_step(data)
b = E_step_old(data)
sum(colSums(a$V == b$V))
nrow(a$V)*ncol(a$V)

data = data0
data = E_step(data)
a = MSU_step(data)
b = MSU_step_old(data)
sum(colSums(a$gamma_y == b$gamma_y))
nrow(a$gamma_y)*ncol(a$gamma_y)


