library(tidyverse)
library(furrr)
plan(multicore, workers = 10)

# Path
bls_path = "./data/BLS/"

# Detect the BLS QCEW files
files_qtrly = list.files(pattern = "q1-q4", path = bls_path, all.files = TRUE)
files_annual = list.files(pattern = "annual", path = bls_path, all.files = TRUE)

# Read the BLS QCEW files
raw_bls_qtrly = files_qtrly %>%
  as.list() %>%
  set_names() %>% 
  future_map(
    ~ read.csv(paste0(bls_path,.), header = TRUE, sep = ",")
  )

raw_bls_annual = files_annual %>%
  as.list() %>%
  set_names() %>% 
  future_map(
    ~ read.csv(paste0(bls_path,.), header = TRUE, sep = ",")
  )

# Save data as .Rda
saveRDS(
  list(
    raw_bls_qtrly = raw_bls_qtrly,
    raw_bls_annual = raw_bls_annual
  ),
  file = "./data/raw_bls_qcew.Rds"
)
