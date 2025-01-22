library(tidyverse)

data = readRDS("./data/qcew_true.Rds")
# data = readRDS("./data/qcew_sup.Rds")
# data = readRDS("./data/qcew_est.Rds")

industry_list = c("3361", "3369", "3371")

output = data %>% 
  filter(industry_code %in% industry_list)

write.csv(output, "./data/qcew_selected.csv")