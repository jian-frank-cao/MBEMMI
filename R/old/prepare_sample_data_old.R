library(tidyverse)
library(furrr)
plan(multicore, workers = 10)
setwd("F:/Caltech/Paul/MBEMMI")

## Read bls and deo data -------------------------------------------------------------
data_bls = readRDS(
  "F:/Caltech/Paul/QCEW/data/data_bls_qcew.Rds"
)

data_deo = readRDS(
  "F:/Caltech/Paul/QCEW/data/data_deo_qcew.Rds"
)

## Merge bls and deo qtrly data ------------------------------------------------------------
# replace null with -999
emplvl_qtrly = data_deo$v2012$qtrly$emplvl_qtrly
emplvl_qtrly[is.na(emplvl_qtrly)] = -999

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
  data_qtrly[3:ncol(data_qtrly)] == -999
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

# replace null with 0
data_qtrly[data_qtrly == -999] = 0

# compute five-digit level
aggregate_6digit = function(data){
  out = case_when(
    is.na(data) ~ sum(data, na.rm = TRUE),
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

# make sure industry_code is consistent
data_qtrly = data_qtrly[match(data_bls$v2012$naics_codes$industry_code, data_qtrly$industry_code),]
null_qtrly = null_qtrly[match(data_bls$v2012$naics_codes$industry_code, null_qtrly$industry_code),]

# get disclosure indicators
disclosure_qtrly = data_bls$v2012$qtrly$disclosure_code


## Merge bls and deo annual data ----------------------------------------------------------
# replace null with -999
emplvl_annual = data_deo$v2012$annual$Average_Monthly_Employment
emplvl_annual[is.na(emplvl_annual)] = -999

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
  data_annual[3:ncol(data_annual)] == -999
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

# replace null with 0
data_qtrly[data_qtrly == -999] = 0


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

# make sure industry_code is consistent
data_annual = data_annual[match(data_bls$v2012$naics_codes$industry_code, data_annual$industry_code),]
null_annual = null_annual[match(data_bls$v2012$naics_codes$industry_code, null_annual$industry_code),]

# locate disclosure indicators
disclosure_annual = data_bls$v2012$annual$disclosure_code


## Picking data sample ------------------------------------------------------------------------
target = c("336", "3361", "3362", "3363", "3364", "3365", "3366", "3369")
ind = (data_qtrly$industry_code %in% target) %>% 
  which
  
sample = list(
  suppressed = list(
    qtrly = data_qtrly[ind,],
    annual = data_annual[ind,]
  ),
  estab = list(
    qtrly = data_bls$v2012$qtrly$qtrly_estabs_count[ind,]
  ),
  true = list(
    qtrly = data_qtrly[ind,],
    annual = data_annual[ind,]
  ),
  disclosure = list(
    qtrly = disclosure_qtrly[ind,],
    annual = disclosure_annual[ind,]
  ),
  null = list(
    qtrly = null_qtrly[ind,],
    annual = null_annual[ind,]
  )
)

sample$suppressed$qtrly[disclosure_qtrly[ind,] == "N" |
                          disclosure_qtrly[ind,] == "-"] = NA
sample$suppressed$annual[disclosure_annual[ind,] == "N" |
                          disclosure_annual[ind,] == "-"] = NA

saveRDS(
  sample,
  "./data/qcew_sample.Rds"
)









  
  
