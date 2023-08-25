library(parallel)
n.cores = detectCores()
library(tidyverse)
library(furrr)
plan(multisession, workers = n.cores - 1)

# Path
deo_path = "./data/DEO/"

# Detect the DEO QCEW files
files_list = list.files(pattern = ".txt", path = deo_path, all.files = TRUE)

# Read the DEO QCEW files
raw_deo = files_list %>%
  as.list() %>%
  set_names() %>% 
  future_map(
    ~ read.csv(
      paste0(deo_path,.),
      header = FALSE,
      sep = ";",
      col.names = c("Year", "Area_Type", "Area_Code", "Ownership_Code",
                    "NAICS_Short_Title", "NAICS_or_BLS_Code", "Confidentiality_Flag",
                    "Units", "Month_1_Employment", "Month_2_Employment",
                    "Month_3_Employment", "Month_4_Employment", "Month_5_Employment",
                    "Month_6_Employment", "Month_7_Employment", "Month_8_Employment",
                    "Month_9_Employment", "Month_10_Employment", "Month_11_Employment",
                    "Month_12_Employment", "Total_Wages", "Taxable_Wages", "Contributions",
                    "Average_Monthly_Employment", "Average_Annual_Wage", "Average_Weekly_Wage",
                    "Month_1_Female_Employment", "Month_2_Female_Employment",
                    "Month_3_Female_Employment", "Month_4_Female_Employment",
                    "Month_5_Female_Employment", "Month_6_Female_Employment",
                    "Month_7_Female_Employment", "Month_8_Female_Employment",
                    "Month_9_Female_Employment", "Month_10_Female_Employment",
                    "Month_11_Female_Employment", "Month_12_Female_Employment", "End")
    )
  )

# Save data as .Rda
saveRDS(
  raw_deo,
  file = "./data/raw_deo_qcew.Rds"
)
