# Parallel Sequential Imputation
PSI_slow <- function(data_qcew, method = "MBEMMI", n_imp = 10, ...){
  # partition data
  blocks <- partition_data(data_qcew)
  
  # generate sobol sequence
  n_time <- nrow(blocks[[3]][[1]][[1]])
  n_block <- sum(sapply(blocks, length))
  sobol_seq <- prepare_sobol(n = n_time * n_block, dim = n_imp * 5)
  list_sobol_imp <- lapply(1:n_imp, function(i){
    start_col <- (i - 1) * 5 + 1
    end_col <- min(i * 5, n_imp * 5)
    sobol_seq[, start_col:end_col, drop = FALSE]
  })
  
  # sequential imputations
  result_imp = NULL
  for (num_imp in 1:n_imp) {
    sobol_imp <- list_sobol_imp[[num_imp]]
    data_imputed <- data_qcew$suppressed$qtrly
    
    for (num_level in 1:length(blocks)) {
      list_data_block <- blocks[[num_level]]
      n_block_level <- length(list_data_block)
      if (n_block_level == 0) next
      
      # prepare sobol sequences for blocks
      list_sobol_block <- lapply(1:n_block_level, function(i){
        start_row <- (i - 1) * n_time + 1
        end_row <- min(i * n_time, n_block_level * n_time)
        sobol_imp[start_row:end_row, , drop = FALSE]
      })
      sobol_imp <- sobol_imp[-(1:n_block_level * n_time), ]
      
      # sequential blocks
      result_level = NULL
      for (num_block in 1:length(list_data_block)) {
        data <- list_data_block[[num_block]]
        sobol_block <- list_sobol_block[[num_block]]
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
        
        # compute missing values if possible
        if (length(code_lower) > 1) {
          n_row_missing = rowSums(is.na(data$detail), na.rm = T)
          if (any(n_row_missing == 1)) {
            ind_row_missing = which(n_row_missing == 1)
            for (ind in ind_row_missing) {
              row = data$detail[ind,]
              row[is.na(row)] = data$q_total[ind,] - sum(row, na.rm = T)
              data$detail[ind,] = row
            }
          }
        }
        
        # impute missing values
        if (length(code_lower) == 1) {
          output = list()
        }else if (!any(n_row_missing > 1)) {
          output = data$detail
        }else{
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
            output = res[["data_final"]][[1]]
          }else if(method == "BMMI"){
            res = safe_MBEMMI(data, n_imputation = 1)
            if (!is.null(res) & length(res) > 0) {
              output = res
            }else{
              output = data$detail
            }
          }else if(method == "EMB"){
            print(paste0("imp:", num_imp, ", code:", code_higher))
            if (num_imp == 8 & code_higher == "33999") {
              res = NULL
            }else if(num_imp == 18 & code_higher == "336"){
              res = NULL
            }else{
              capture.output(res <- safe_EMB(data, n_imputation = 1),
                             file = "/dev/null")
            }
            # res <- safe_EMB(data, n_imputation = 1)
            if (!is.null(res) & length(res) > 0) {
              output = res
            }else{
              output = data$detail
            }
          }
        }
        result_level[[names(list_data_block)[num_block]]] = output
        
      }
      
      # Loop over each name in data_level
      for (item in names(list_data_block)) {
        
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
    result_imp[[paste0("imputation-", num_imp)]] = data_imputed
  }
  return(result_imp)
}


safe_MBEMMI <- function(data, n_imputation) {
  tryCatch(
    {
      output = BMMI(data, n_imputation)[[1]]
      return(output)
    },
    error = function(e) {
      message("Error: ", e$message)
      return(NULL) 
    }
  )
}


safe_EMB <- function(data, n_imputation) {
  tryCatch(
    {
      output = R.utils::withTimeout(
        EMB(data, n_imputation)[[1]],
        timeout = 10,  # Timeout in seconds
        onTimeout = "error"
      )
      return(output)
    },
    error = function(e) {
      message("Error: ", e$message)
      return(NULL) 
    }
  )
}
