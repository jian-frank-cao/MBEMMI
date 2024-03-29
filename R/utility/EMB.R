## Packages --------------------------------------------------------------------
# library(Amelia)

## EMB -------------------------------------------------------------------------
EMB = function(data_raw, n_imputation = 10,
               polytime = 1){
  # prepare data
  data = cbind(data.frame(ts = 1:nrow(data_raw$detail)),
               data_raw$detail)
  
  result <- tryCatch(
    {
      res <- Amelia::amelia(data, m = n_imputation,
                            ts = 'ts', polytime = polytime)
      res = res[["imputations"]]
      for (i in 1:length(res)) {
        res[[i]] = res[[i]][-1]
      }
      res
    },
    error = function(e) NULL
  )

  return(result)
}
