## Packages --------------------------------------------------------------------
# sudo spctl --global-disable
library(parallel)
n.cores = detectCores()
library(tidyverse)
library(furrr)
plan(multisession, workers = n.cores - 2)
set.seed(1234)

source("./R/utility/MBEMMI.R")
source("./R/utility/BMMI.R")
source("./R/utility/EMB.R")

## Read data -------------------------------------------------------------------
data_qcew = readRDS("./data/qcew_2digit_fixed.Rds")

## Functions -------------------------------------------------------------------
# Functions to partition data
partition_data = function(data_qcew){
  
  blocks = 2:5 %>% 
    as.list %>% 
    set_names(paste0("Digit-", 2:5)) %>% 
    future_map(
      ~{
        digit = .
        data = data_qcew$suppressed$qtrly %>% 
          filter(nchar(industry_code) %in% c(digit, digit + 1)) %>% 
          mutate(
            group = substr(industry_code, 1, digit),
            n_digit = nchar(industry_code)
          )
        data$n_missing = is.na(data) %>% rowSums(.)
        groups = data$group[data$n_missing > 0 &
                              data$n_digit == digit + 1] %>% unique
        data = data %>% 
          filter(group %in% groups) %>% 
          dplyr::select(-c("n_missing", "n_digit")) %>% 
          group_by(group) %>%
          group_split(
            .,
            .keep = FALSE
          ) %>% 
          set_names(groups)
        a = data %>% 
          future_map(
            ~{
              block = .
              ind_digit = nchar(block$industry_code) == digit
              ind_digit_1 = nchar(block$industry_code) == digit + 1
              code_higher = block$industry_code[ind_digit]
              code_lower = block$industry_code[ind_digit_1]
              list(
                detail = block %>%
                  filter(industry_code %in% code_lower) %>% 
                  dplyr::select(-c("industry_code", "industry_title")) %>% 
                  t %>% 
                  data.frame %>% 
                  `colnames<-`(code_lower),
                auxiliary = data_qcew$estab$qtrly %>%
                  filter(industry_code %in% code_lower) %>% 
                  dplyr::select(-c("industry_code", "industry_title")) %>% 
                  t %>% 
                  data.frame %>% 
                  `colnames<-`(code_lower),
                q_total = block %>%
                  filter(industry_code %in% code_higher) %>% 
                  dplyr::select(-c("industry_code", "industry_title")) %>% 
                  t %>% 
                  data.frame %>% 
                  `colnames<-`(code_higher),
                a_total = data_qcew$suppressed$annual %>%
                  filter(industry_code %in% code_lower) %>% 
                  dplyr::select(-c("industry_code", "industry_title")) %>% 
                  # mutate_all(~ .*4) %>% 
                  t %>% 
                  data.frame %>% 
                  `colnames<-`(code_lower),
                a_q_total = data_qcew$suppressed$annual %>%
                  filter(industry_code %in% code_higher) %>% 
                  dplyr::select(-c("industry_code", "industry_title")) %>% 
                  # mutate_all(~ .*4) %>% 
                  t %>% 
                  data.frame %>% 
                  `colnames<-`(code_higher)
              )
            }
          )
      }
    )
  
  return(blocks)
}

# Parallel Sequential Imputation
PSI <- function(data_qcew, method = "MBEMMI", n_imp = 10,
                parallel_imp = TRUE, parallel_level = FALSE, ...){
  # define map type
  fun_map_imp <- ifelse(parallel_imp, future_map, map)
  fun_map_level <- ifelse(parallel_level, future_map2, map2)
  
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
  
  # parallel imputations
  result_imp <- list_sobol_imp %>% 
    set_names(paste0("imputation-", 1:n_imp)) %>% 
    fun_map_imp(
      .options=furrr_options(seed = TRUE),
      ~{
        sobol_imp <- .
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
          result_level <- fun_map_level(
            list_data_block,
            list_sobol_block,
            ~{
              data <- .x
              sobol_block <- .y
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
                  res[["data_final"]][[1]]
                }else if(method == "BMMI"){
                  BMMI(data, n_imputation = 1)[[1]]
                }else if(method == "EMB"){
                  res <- EMB(data, n_imputation = 1)[[1]]
                  if (length(res) == 1) {
                    res
                  }else{
                    data$detail
                  }
                }
              }else{
                list()
              }
            }
          )
          
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
        data_imputed
      }
    )
  
  return(result_imp)
}


## Clean missing > 60% ---------------------------------------------------------
blocks = partition_data(data_qcew)

industry_60 = NULL
for (level in blocks) {
  for (block in level) {
    n_missing = colSums(is.na(block$detail))
    if (sum(n_missing >=12) > 0) {
      industry_60 = c(industry_60, colnames(block$detail))
    }
  }
}

ind_industry_60 = which(data_qcew$suppressed$qtrly$industry_code %in%
                          industry_60)
sum(data_qcew$suppressed$qtrly$industry_code[ind_industry_60] != 
      data_qcew$true$qtrly$industry_code[ind_industry_60])
  
data_qcew$suppressed$qtrly[ind_industry_60, ] = 
  data_qcew$true$qtrly[ind_industry_60, ]

## Compute missing values if possible ------------------------------------------
blocks = partition_data(data_qcew)

for (level in blocks) {
  for (block in level) {
    code_higher <- colnames(block$q_total)
    code_lower <- colnames(block$detail)
    
    if (any(is.na(block$q_total))) {  
      block$q_total <- data_qcew$suppressed$qtrly %>%
        filter(industry_code %in% code_higher) %>%
        dplyr::select(-c("industry_code", "industry_title")) %>%
        t() %>%
        as.data.frame() %>%
        dplyr::rename_all(~code_higher)
    }
    
    if (length(code_lower) > 1) {
      n_row_missing = rowSums(is.na(block$detail), na.rm = T)
      if (any(n_row_missing == 1)) {
        ind_row_missing = which(n_row_missing == 1)
        for (ind in ind_row_missing) {
          row = block$detail[ind,]
          industry = colnames(row)[which(is.na(row))]
          ind_industry = which(data_qcew$suppressed$qtrly$industry_code == industry)
          data_qcew$suppressed$qtrly[ind_industry, ind + 2] =
            block$q_total[ind,] - sum(row, na.rm = T)
        }
      }
    }
  }
}


## Run -------------------------------------------------------------------------
source("./R/utility/PSI.R")

set.seed(1234)
result_MBEMMI = PSI_slow(data_qcew, method = "MBEMMI", n_imp = 10)
saveRDS(result_MBEMMI, paste0("./data/imputations_MBEMMI_florida.Rds"))

set.seed(1234)
result_BMMI = PSI_slow(data_qcew, method = "BMMI", n_imp = 10)
saveRDS(result_BMMI, paste0("./data/imputations_BMMI_florida.Rds"))

set.seed(1234)
result_EMB = PSI_slow(data_qcew, method = "EMB", n_imp = 30)
select_EMB = NULL
for (i in 1:length(result_EMB)) {
  if (sum(is.na(result_EMB[[i]])) == 50) {
    select_EMB = c(select_EMB, i)
  }
}
saveRDS(result_EMB[select_EMB[1:10]],
        paste0("./data/imputations_EMB_florida.Rds"))

