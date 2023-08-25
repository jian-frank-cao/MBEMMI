## Packages --------------------------------------------------------------------
# library(Amelia)

## EMB -------------------------------------------------------------------------
EMB = function(data_raw, n_imputation = 10,
               polytime = 1){
  # prepare data
  data = cbind(data.frame(ts = 1:nrow(data_raw$detail)),
               data_raw$detail)
  
  result = Amelia::amelia(data[-6],
                          m = n_imputation,
                          ts = 'ts',
                          polytime = polytime)[["imputations"]]
  
  for (i in 1:length(result)) {
    result[[i]] = result[[i]][-1]
  }
  
  return(result)
}