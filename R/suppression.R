## Packages ------------------------------------------------------------------------------------
library(tidyverse)
library(furrr)
plan(multiprocess, workers = 10)
set.seed(1234)

## Read data -----------------------------------------------------------------------------------
data_qcew = readRDS("./data/qcew_2digit_fixed.Rds")


## Random Primary Suppression ------------------------------------------------------------------
rnd_primary_sup = function(data, digit = 6,
                         value_sup_rate = 0.05, 
                         row_sup_rate = 0.20){
  # target level
  data_target = data %>%
    filter(nchar(industry_code) == digit)
  
  # rank the industry by means
  ranking = data.frame(industry_code = data_target$industry_code,
                       mean = rowMeans(data_target[-c(1,2)])) %>% 
    arrange(mean)
  
  # randomly select rows using geometric distribution (prob = 0.0068)
  n_row = nrow(data_target)
  n_row_sup = as.integer(n_row * row_sup_rate)
  geom_prob = dgeom(1:n_row, 0.0068)
  ranking_sample = sample(n_row, n_row_sup, replace = FALSE, prob = geom_prob)
  industry_sample = ranking$industry_code[ranking_sample]
  data_sample = data_target %>% filter(industry_code %in% industry_sample)
  
  # random suppression
  n_sup = nrow(data_target) * (ncol(data_target) - 2) * value_sup_rate
  n_y_sample = nrow(data_sample) * (ncol(data_sample) - 2)
  ind_sup = sample(1:n_y_sample, n_sup, replace = FALSE)
  vector_sup = rep(FALSE, n_y_sample)
  vector_sup[ind_sup] = TRUE
  vector_sup = c(rep(FALSE, nrow(data_sample) * 2), vector_sup)
  matrix_sup = matrix(vector_sup, nrow = nrow(data_sample))
  data_sample[matrix_sup] = NA_real_
  
  # map to data
  data[match(data_sample$industry_code, data$industry_code), 3:ncol(data)] = 
    data_sample[, 3:ncol(data_sample)]
  
  return(data)
}

## Secondary Suppresion ------------------------------------------------------------------------
split_blocks = function(data, data_annual, level){
  data_list = data %>% 
    filter(nchar(industry_code) %in% level:(level+1)) %>% 
    mutate(group = substr(industry_code, 1, level)) %>% 
    group_by(group) %>% 
    group_split(.keep = FALSE)
  
  data_annual_list = data_annual %>% 
    filter(nchar(industry_code) %in% level:(level+1)) %>% 
    mutate(group = substr(industry_code, 1, level)) %>% 
    group_by(group) %>% 
    group_split(.keep = FALSE)
  
  res = future_map2(
    data_list,
    data_annual_list,
    ~{
      qtrly = .x
      annual = .y
      list(qtrly_total = qtrly %>%
             filter(nchar(industry_code) == level),
           qtrly_detail = qtrly %>%
             filter(nchar(industry_code) == level + 1),
           annual_total = annual %>%
             filter(nchar(industry_code) == level),
           annual_detail = annual %>%
             filter(nchar(industry_code) == level + 1)
           )
    }
  )
  
  return(res)
}

sup_block = function(data){
  qtrly_total = data$qtrly_total
  qtrly_detail = data$qtrly_detail
  annual_total = data$annual_total
  annual_detail = data$annual_detail
  
  n = nrow(qtrly_detail)
  ind_sup = is.na(qtrly_detail)
  
  if (sum(ind_sup) == 0) {
    return(data)
  }
  
  # NAICS code dimension
  if (n == 1) {
    qtrly_total[ind_sup] = NA
  }else{
    ind_cols = which(colSums(ind_sup) == 1)
    for (col in ind_cols) {
      col_sup = ind_sup[,col]
      ind_rows = which(col_sup == FALSE)
      if (length(ind_rows) == 1) {
        qtrly_detail[ind_rows, col] = NA
      }else{
        year = (col+1) %/% 4
        block = qtrly_detail[ind_rows, (year*4-1):(year*4+2)]
        block_sup = is.na(block)
        exist_sup = which(rowSums(block_sup) > 0)
        if (length(exist_sup) > 0) {
          ind_rows = ind_rows[exist_sup]
        }
        temp_col = qtrly_detail[ind_rows, col]
        ind = which(temp_col == min(temp_col))
        qtrly_detail[ind_rows[ind], col] = NA
      }
    }
  }
  
  # time dimension
  for (year in 1:5) {
    block = qtrly_detail[(year*4-1):(year*4+2)]
    block_sup = is.na(block)
    ind_rows = which(rowSums(block_sup) == 1)
    annual_detail[ind_rows, year + 2] = NA
  }
  
  data$qtrly_total = qtrly_total
  data$qtrly_detail = qtrly_detail
  data$annual_total = annual_total
  data$annual_detail = annual_detail
  
  return(data)
}

secondary_sup = function(data_qtrly, data_annual){
  finished = FALSE
  count = 1
  
  while (!finished) {
    print(count)
    sup_qtrly = is.na(data_qtrly)
    sup_annual = is.na(data_annual)
    
    for (level in 5:2) {
      data_list = split_blocks(data_qtrly, data_annual, level)
      data_list = data_list %>% 
        future_map(
          ~{
            sup_block(.)
          }
        )
      qtrly_total = lapply(data_list, "[[", "qtrly_total") %>% 
        do.call("rbind", .)
      qtrly_detail = lapply(data_list, "[[", "qtrly_detail") %>% 
        do.call("rbind", .)
      annual_total = lapply(data_list, "[[", "annual_total") %>% 
        do.call("rbind", .)
      annual_detail = lapply(data_list, "[[", "annual_detail") %>% 
        do.call("rbind", .)
      
      data_qtrly[match(qtrly_total$industry_code, data_qtrly$industry_code),] = 
        qtrly_total
      data_qtrly[match(qtrly_detail$industry_code, data_qtrly$industry_code),] = 
        qtrly_detail
      data_annual[match(annual_total$industry_code, data_annual$industry_code),] = 
        annual_total
      data_annual[match(annual_detail$industry_code, data_annual$industry_code),] = 
        annual_detail
    }
    
    qtrly_finished = sum(sup_qtrly != is.na(data_qtrly)) == 0
    annual_finished = sum(sup_annual != is.na(data_annual)) == 0
    if (qtrly_finished & annual_finished) {
      finished = TRUE
    }
    count = count + 1
  }
  
  print("Finished!")
  return(list(qtrly = data_qtrly,
              annual = data_annual))
}

## Generate 50 suppressed data sets ------------------------------------------------------------
set.seed(1234)
n_data_set = 50
seed_list = sample(10000, n_data_set, replace = FALSE)

data_list = NULL
for (i in 1:n_data_set) {
  print(paste0("Random suppression ",i))
  data_primary = rnd_primary_sup(data_qcew$true$qtrly)
  data_list[[i]] = data_qcew
  data_list[[i]][["suppressed"]] = secondary_sup(data_primary,
                                                 data_list[[i]]$true$annual)
}

saveRDS(data_list, "./data/qcew_rnd_suppression.Rds")

