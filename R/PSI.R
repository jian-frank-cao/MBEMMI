## Packages ------------------------------------------------------------------------------------
# sudo spctl --global-disable
library(checkpoint)
checkpoint("2020-09-15")

library(tidyverse)
library(furrr)
plan(multiprocess, workers = 10)
set.seed(1234)

source("./R/MBEMMI_20210929.R")
source("./R/BMMI.R")
source("./R/EMB.R")

## Read data -----------------------------------------------------------------------------------
# data_qcew = readRDS("./data/qcew_2digit_fixed.Rds")
data_qcew = readRDS("./data/qcew_rnd_suppression.Rds")[[1]]

## Functions to partition data -----------------------------------------------------------------
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
        data %>% 
          future_map(
            ~{
              block = .
              code_higher = block$industry_code[nchar(block$industry_code) == digit]
              code_lower = block$industry_code[nchar(block$industry_code) == digit + 1]
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

## PSI ------------------------------------------------------------------------------------------
PSI = function(data_qcew, method = "MBEMMI",
               n_imputation = 10, ...){
  
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
  
  return(result)
}


## Run -----------------------------------------------------------------------------------------
result = PSI(data_qcew)

for (i in 1:length(blocks$`Digit-3`)) {
  print(i)
  data = blocks$`Digit-3`[[i]]
  a = MBEMMI(data, n_imputation = 2)
}


method = "MBEMMI"
imputation = NULL
n_imputation = 10

blocks = partition_data(data_qcew)


imputation = 1:n_imputation %>% 
  as.list %>% 
  set_names(paste0("imputation-", 1:n_imputation)) %>% 
  future_map(
    ~{
      data_imputed = data_qcew$suppressed$qtrly
      
      for (level in 1:length(blocks)) {
        data_level = blocks[[level]]
        result = NULL
        
        for (item in names(data_level)[1:length(data_level)]) {
          if (item != "22111") {
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
                converged = FALSE
                while(converged == FALSE){
                  res = NULL
                  tryCatch({res = MBEMMI(data, n_imputation = 1)},
                           error = function(err) {
                             print(paste("MY ERROR: ", err))
                           })
                  if (!is.null(res)) {
                    if (res[[1]]$converged == TRUE) {
                      result[[item]] = res[[1]][["data_final"]][[1]]
                      converged = TRUE # solve sobol seq problem + add qmc at simulation
                    }
                  }
                }
              }else if(method == "BMMI"){
                BMMI(data, n_imputation = 1)[[1]]
              }else if(method == "EMB"){
                EMB(data, n_imputation = 1)[[1]]
              }
            }else{
              list()
            }
          }
        }
        
        
        for (item in names(data_level)) {
          code_lower = colnames(data_level[[item]]$detail)
          code_higher = colnames(data_level[[item]]$q_total)
          if (item != "22111") {
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
      rows = which(substr(data_imputed$industry_code, 1, 5) == "22111")
      data_imputed[rows, ] = data_qcew$true$qtrly[rows, ]
      data_imputed
    })
  
saveRDS(imputation, "./data/QCEW_imputations_20211006.Rds")




# # which(names(data_level) == item)
# # 4-32, 
# 
# for (item in names(data_level)[1:length(data_level)]) {
#   if (item != "22111") {
#     data = data_level[[item]]
#     code_higher = colnames(data$q_total)
#     code_lower = colnames(data$detail)
#     if (any(is.na(data$q_total))) {  
#       data$q_total = data_imputed %>%
#         filter(industry_code %in% code_higher) %>% 
#         dplyr::select(-c("industry_code", "industry_title")) %>% 
#         t %>% 
#         data.frame %>% 
#         `colnames<-`(code_higher)
#     }
#     if (length(code_lower) > 1) {
#       if (method == "MBEMMI") {
#         converged = FALSE
#         while(converged == FALSE){
#           res = NULL
#           tryCatch({res = MBEMMI(data, n_imputation = 1)},
#                    error = function(err) {
#                      print(paste("MY ERROR: ", err))
#                    })
#           if (!is.null(res)) {
#             if (res[[1]]$converged == TRUE) {
#               result[[item]] = res[[1]][["data_final"]][[1]]
#               converged = TRUE # solve sobol seq problem + add qmc at simulation
#             }
#           }
#         }
#       }else if(method == "BMMI"){
#         BMMI(data, n_imputation = 1)[[1]]
#       }else if(method == "EMB"){
#         EMB(data, n_imputation = 1)[[1]]
#       }
#     }else{
#       list()
#     }
#   }
# }
# 
# 
# for (item in names(data_level)) {
#   code_lower = colnames(data_level[[item]]$detail)
#   code_higher = colnames(data_level[[item]]$q_total)
#   if (item != "22111") {
#     if (length(code_lower) > 1) {
#       data_imputed[data_imputed$industry_code %in% code_lower,
#                    3:ncol(data_imputed)] = result[[item]] %>% t
#     }else{
#       data_imputed[data_imputed$industry_code %in% code_lower,
#                    3:ncol(data_imputed)] = data_imputed[data_imputed$industry_code %in% code_higher,
#                                                         3:ncol(data_imputed)]
#     }
#   }
# }
# 
# rows = which(substr(data_imputed$industry_code, 1, 5) == "22111")
# data_imputed[rows, ] = data_qcew$true$qtrly[rows, ]
# 
# 
# 
# 
# 
# 
# m = 2
# imputation[[m]] = data_imputed
# 





