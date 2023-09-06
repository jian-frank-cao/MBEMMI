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
