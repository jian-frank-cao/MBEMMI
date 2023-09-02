library(parallel)
n.cores = detectCores()
library(tidyverse)
library(furrr)
plan(multisession, workers = n.cores - 1)

## Read bls and deo data -------------------------------------------------------
data_bls = readRDS(
  "./data/data_bls_qcew.Rds"
)

data_deo = readRDS(
  "./data/data_deo_qcew.Rds"
)

## Merge bls and deo qtrly data ------------------------------------------------
# replace null with -Inf
emplvl_qtrly = data_deo$v2012$qtrly$emplvl_qtrly
emplvl_qtrly[is.na(emplvl_qtrly)] = -Inf

# merge
data_qtrly = merge(
  data_bls$v2012$naics_codes,
  emplvl_qtrly,
  by.x = "industry_code",
  by.y = "NAICS_or_BLS_Code",
  all = TRUE
)

# remove redundant industry_codes from deo
data_qtrly = data_qtrly %>% 
  filter(!is.na(industry_title)) %>% 
  select(-NAICS_Short_Title)

# get null indicators
determine_null = function(data){
  out = case_when(
    (is.na(data) & (sum(data, na.rm = TRUE) == (length(data) - 1))) ~ TRUE,
    (is.na(data) & (sum(data, na.rm = TRUE) < (length(data) - 1))) ~ FALSE,
    TRUE ~ data
  )
}

null_qtrly = cbind(
  data_qtrly[1:2],
  data_qtrly[3:ncol(data_qtrly)] == -Inf
)

null_qtrly = null_qtrly %>% 
  mutate(
    group_5digit = case_when(
      nchar(industry_code) > 4 ~ industry_code %>% substr(., 1, 5),
      TRUE ~ "NA"
    )
  ) %>% 
  group_by(group_5digit) %>% 
  mutate_at(
    vars(colnames(null_qtrly)[3:22]),
    ~ determine_null(.)
  ) %>% 
  ungroup %>% 
  select(-c("group_5digit"))

# compute five-digit level
aggregate_6digit = function(data){
  out = case_when(
    is.na(data)  ~ sum(data, na.rm = TRUE),
    TRUE ~ data
  )
  return(out)
}

data_qtrly = data_qtrly %>% 
  mutate(
    group_5digit = case_when(
      nchar(industry_code) > 4 ~ industry_code %>% substr(., 1, 5),
      TRUE ~ "NA"
    )
  ) %>% 
  group_by(group_5digit) %>% 
  mutate_at(
    vars(colnames(data_qtrly)[3:22]),
    ~ aggregate_6digit(.)
  ) %>% 
  ungroup %>% 
  select(-c("group_5digit"))

# compute null
data_qtrly[data_qtrly == -Inf] = NA_real_
data_qtrly[data_qtrly$industry_code == "31123", c(3,4,5,6)] = 0

for (col in colnames(data_qtrly)[-c(1,2)]) {
  null_ind = which(is.na(data_qtrly[col]))
  for (ind in null_ind) {
    code = data_qtrly$industry_code[ind]
    ndigit = nchar(code)
    higher_code = substr(code, 1, ndigit - 1)
    block_ind = which((substr(data_qtrly$industry_code, 1,
                              ndigit - 1) == higher_code) &
                        (nchar(data_qtrly$industry_code) < ndigit + 1))
    total = data_qtrly[block_ind[1], col] %>% unlist
    detail = data_qtrly[block_ind[-1], col] %>% unlist
    if (sum(is.na(detail)) == 1 & !is.na(total)) {
      data_qtrly[ind, col] = round(total - sum(detail, na.rm = TRUE), 10)
    }
  }
}

# make sure industry_code is consistent
data_qtrly = data_qtrly[match(data_bls$v2012$naics_codes$industry_code,
                              data_qtrly$industry_code),]
null_qtrly = null_qtrly[match(data_bls$v2012$naics_codes$industry_code,
                              null_qtrly$industry_code),]

# get disclosure indicators
disclosure_qtrly = data_bls$v2012$qtrly$disclosure_code

# round to integer
data_qtrly[-c(1,2)] = round(data_qtrly[-c(1,2)], 0)


## Merge bls and deo annual data -----------------------------------------------
# replace null with -Inf
emplvl_annual = data_deo$v2012$annual$Average_Monthly_Employment
emplvl_annual[is.na(emplvl_annual)] = -Inf

# merge
data_annual = merge(
  data_bls$v2012$naics_codes,
  emplvl_annual,
  by.x = "industry_code",
  by.y = "NAICS_or_BLS_Code",
  all = TRUE
)

# remove redundant industry_codes from deo
data_annual = data_annual %>% 
  filter(!is.na(industry_title)) %>% 
  select(-NAICS_Short_Title)

# get null indicators
null_annual = cbind(
  data_annual[1:2],
  data_annual[3:ncol(data_annual)] == -Inf
)

null_annual = null_annual %>% 
  mutate(
    group_5digit = case_when(
      nchar(industry_code) > 4 ~ industry_code %>% substr(., 1, 5),
      TRUE ~ "NA"
    )
  ) %>% 
  group_by(group_5digit) %>% 
  mutate_at(
    vars(colnames(null_annual)[3:7]),
    ~ determine_null(.)
  ) %>% 
  ungroup %>% 
  select(-c("group_5digit"))

# compute five-digit level
data_annual = data_annual %>% 
  mutate(
    group_5digit = case_when(
      nchar(industry_code) > 4 ~ industry_code %>% substr(., 1, 5),
      TRUE ~ "NA"
    )
  ) %>% 
  group_by(group_5digit) %>% 
  mutate_at(
    vars(colnames(data_annual)[3:7]),
    ~ aggregate_6digit(.)
  ) %>% 
  ungroup %>% 
  select(-group_5digit)

# compute null
data_annual[data_annual == -Inf] = NA_real_
data_annual[data_annual$industry_code == "31123", 3] = 0

for (col in colnames(data_annual)[-c(1,2)]) {
  null_ind = which(is.na(data_annual[col]))
  for (ind in null_ind) {
    code = data_annual$industry_code[ind]
    ndigit = nchar(code)
    higher_code = substr(code, 1, ndigit - 1)
    block_ind = which((substr(data_annual$industry_code, 1,
                              ndigit - 1) == higher_code) &
                        (nchar(data_annual$industry_code) < ndigit + 1))
    total = data_annual[block_ind[1], col] %>% unlist
    detail = data_annual[block_ind[-1], col] %>% unlist
    if (sum(is.na(detail)) == 1 & !is.na(total)) {
      data_annual[ind, col] = round(total - sum(detail, na.rm = TRUE), 10)
    }
  }
}

data_annual[data_annual<0] = 0

# make sure industry_code is consistent
data_annual = data_annual[match(data_bls$v2012$naics_codes$industry_code,
                                data_annual$industry_code),]
null_annual = null_annual[match(data_bls$v2012$naics_codes$industry_code,
                                null_annual$industry_code),]

# locate disclosure indicators
disclosure_annual = data_bls$v2012$annual$disclosure_code

data_annual[-c(1,2)] = data_annual[-c(1,2)] *4


## Create alternative data_annual ----------------------------------------------
sum_operator = kronecker(diag(1, 5), rep(1,4))
data_annual_alt = cbind(data_qtrly[c(1,2)],
                        as.matrix(data_qtrly[-c(1,2)]) %*% sum_operator)
colnames(data_annual_alt) = colnames(data_annual)


## Save processed data ---------------------------------------------------------
data = list(
  suppressed = list(
    qtrly = data_qtrly,
    annual = data_annual_alt,
    annual_orig = data_annual
  ),
  estab = list(
    qtrly = data_bls$v2012$qtrly$qtrly_estabs_count
  ),
  true = list(
    qtrly = data_qtrly,
    annual = data_annual_alt,
    annual_orig = data_annual
  ),
  disclosure = list(
    qtrly = disclosure_qtrly,
    annual = disclosure_annual
  ),
  null = list(
    qtrly = null_qtrly,
    annual = null_annual
  )
)

data$suppressed$qtrly[disclosure_qtrly == "N" |
                        disclosure_qtrly == "-"] = NA
data$suppressed$annual[disclosure_annual == "N" |
                         disclosure_annual == "-"] = NA
data$suppressed$annual_orig[disclosure_annual == "N" |
                              disclosure_annual == "-"] = NA

# saveRDS(
#   data,
#   "./data/qcew_ready.Rds"
# )


## Handle 2-digit industries ---------------------------------------------------
target_codes = c("31","32","33","44","45","48","49")
target_titles = c("Manufacturing 01","Manufacturing 02","Manufacturing 03",
                  "Retail trade 01","Retail trade 02",
                  "Transportation and Warehousing 01",
                  "Transportation and Warehousing 02")

handle_2_digit = function(data){
  for (i in 1:length(target_codes)) {
    code = target_codes[i]
    title = target_titles[i]
    industry = data[1,]
    industry$industry_code = code
    industry$industry_title = title
    target = data %>%
      filter(nchar(industry_code) == 3 &
               substr(industry_code, 1, 2) == code)
    industry[3:ncol(data)] = t(colSums(target[3:ncol(data)]))
    data = rbind(data, industry)
  }
  
  data = data %>%
    filter(!grepl("-", industry_code, fixed = TRUE)) %>% 
    arrange(industry_code)
  
  return(data)
}

handle_2_digit_v2 = function(data){
  temp = data[1:7,]
  temp$industry_code = target_codes
  temp$industry_title = target_titles
  data = rbind(data, temp) %>% 
    filter(!grepl("-", industry_code, fixed = TRUE)) %>% 
    arrange(industry_code)
  
  return(data)
}

data$estab$qtrly = handle_2_digit(data$estab$qtrly)
data$true$qtrly = handle_2_digit(data$true$qtrly)
data$true$annual = handle_2_digit(data$true$annual)
data$true$annual_orig = handle_2_digit(data$true$annual_orig)
data$disclosure$qtrly = handle_2_digit_v2(data$disclosure$qtrly)
data$disclosure$annual = handle_2_digit_v2(data$disclosure$annual)
data$null$qtrly = handle_2_digit_v2(data$null$qtrly)
data$null$annual = handle_2_digit_v2(data$null$annual)

data$suppressed = data$true
data$suppressed$qtrly[data$disclosure$qtrly == "N" |
                        data$disclosure$qtrly == "-"] = NA
data$suppressed$annual[data$disclosure$annual == "N" |
                         data$disclosure$annual == "-"] = NA
data$suppressed$annual_orig[data$disclosure$annual == "N" |
                              data$disclosure$annual == "-"] = NA

saveRDS(
  data,
  "./data/qcew_2digit_fixed.Rds"
)


# ## Picking data sample -------------------------------------------------------
# target = c("336", "3361", "3362", "3363", "3364", "3365", "3366", "3369")
# ind = (data_qtrly$industry_code %in% target) %>% 
#   which
#   
# sample = list(
#   suppressed = list(
#     qtrly = data_qtrly[ind,],
#     annual = data_annual[ind,]
#   ),
#   estab = list(
#     qtrly = data_bls$v2012$qtrly$qtrly_estabs_count[ind,]
#   ),
#   true = list(
#     qtrly = data_qtrly[ind,],
#     annual = data_annual[ind,]
#   ),
#   disclosure = list(
#     qtrly = disclosure_qtrly[ind,],
#     annual = disclosure_annual[ind,]
#   ),
#   null = list(
#     qtrly = null_qtrly[ind,],
#     annual = null_annual[ind,]
#   )
# )
# 
# sample$suppressed$qtrly[disclosure_qtrly[ind,] == "N" |
#                           disclosure_qtrly[ind,] == "-"] = NA
# sample$suppressed$annual[disclosure_annual[ind,] == "N" |
#                           disclosure_annual[ind,] == "-"] = NA
# 
# saveRDS(
#   sample,
#   "./data/qcew_sample.Rds"
# )
# 








  
  
