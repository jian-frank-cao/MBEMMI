# data = t(data_qcew$suppressed$qtrly)
# df_industry = data[,1:2]
# data = data.frame(data[,-c(1,2)])
data = data_qcew$suppressed$qtrly


# Function to calculate the longest consecutive NA streak in a vector
longest_na_streak <- function(x) {
  rle(is.na(x))$lengths[which.max(rle(is.na(x))$values)] %||% 0
}

# Apply the function to each row of the dataframe
library(purrr) # To use %||% operator
consecutive_missing <- apply(data, 1, longest_na_streak)
consecutive_missing[consecutive_missing == 22] = 0

total_missing = rowSums(is.na(data))

selection1 = which(total_missing <= 8 & total_missing != 0)

industry1 = data[selection1, ]

# industry1 = industry1 %>% 
#   mutate(total_missing = total_missing[selection1],
#          consecutive_missing = consecutive_missing[selection1])

estab1 = data_qcew[["estab"]][["qtrly"]][selection1,]

true1 = data_qcew[["true"]][["qtrly"]][selection1,]

write.csv(industry1, "./data/QCEW_less_than_40_suppressed.csv")
write.csv(estab1, "./data/QCEW_less_than_40_establishment.csv")
write.csv(true1, "./data/QCEW_less_than_40_true.csv")


## two digit level and 333242 --------------------------------------------------
data = data_qcew$true$qtrly %>%
  filter(nchar(industry_code) == 2 | industry_code == "333242")

Manufacturing = data %>%
  filter(industry_code %in% c("31", "32", "33")) %>% 
  mutate(industry_code = 0, industry_title = 0) %>% 
  colSums(.)

Retail = data %>%
  filter(industry_code %in% c("44", "45")) %>% 
  mutate(industry_code = 0, industry_title = 0) %>% 
  colSums(.)

Transportation = data %>%
  filter(industry_code %in% c("48", "49")) %>% 
  mutate(industry_code = 0, industry_title = 0) %>% 
  colSums(.)

out = data %>% 
  filter(!(industry_code %in% c("10", "31", "32", "33",
                                "44", "45", "48", "49", "99"))) %>% 
  rbind(Manufacturing, Retail, Transportation, .)

out$industry_code[1:3] = c("31-33 (exclude 333242)", "44-45", "48-49")
out$industry_title[1:3] = c("Manufacturing (exclude 333242)", "Retail",
                            "Transportation")

out[1,3:22] = out[1,3:22] - out[8,3:22]

write.csv(out, "./data/GVAR_2digit.csv")


## two digit level and 3361 --------------------------------------------------
data = data_qcew$true$qtrly %>%
  filter(nchar(industry_code) == 2 | industry_code == "3361")

Manufacturing = data %>%
  filter(industry_code %in% c("31", "32", "33")) %>% 
  mutate(industry_code = 0, industry_title = 0) %>% 
  colSums(.)

Retail = data %>%
  filter(industry_code %in% c("44", "45")) %>% 
  mutate(industry_code = 0, industry_title = 0) %>% 
  colSums(.)

Transportation = data %>%
  filter(industry_code %in% c("48", "49")) %>% 
  mutate(industry_code = 0, industry_title = 0) %>% 
  colSums(.)

out = data %>% 
  filter(!(industry_code %in% c("10", "31", "32", "33",
                                "44", "45", "48", "49", "99"))) %>% 
  rbind(Manufacturing, Retail, Transportation, .)

out$industry_code[1:3] = c("31-33 (exclude 3361)", "44-45", "48-49")
out$industry_title[1:3] = c("Manufacturing (exclude 3361)", "Retail",
                            "Transportation")

out[1,3:22] = out[1,3:22] - out[8,3:22]

write.csv(out, "./data/GVAR_2digit_v2.csv")
