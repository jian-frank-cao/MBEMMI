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
