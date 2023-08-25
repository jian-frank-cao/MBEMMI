## Packages --------------------------------------------------------------------
# library(MCMCpack)
# library(Matrix)
# library(bayesm)
# library(MASS)
# library(mvtnorm)

## Function to genererate from full conditional of missing data. ---------------
my.missing <- function(z.star,x.star,V,missing.indeces,H,N)
{
  ind <- which(is.na(z.star[(4*N+5):(5*N+4)]))                          
  if (length(ind)>0) z.star <- z.star[-(ind+4*(N+1))]
  
  # The conditional mean of z given mu, sigma2.i, and W.i is:
  mu <- H %*% x.star
  
  # The covariance matrix of z conditional on mu, sigma2.i, and W.i is:
  Sigma <- H %*% V %*% t(H)
  
  # The vector of indices of the observed data is:
  observed.indeces <- (1:nrow(H))[-missing.indeces]
  
  # The partitions of the Sigma matrix according to missing and observed data are:
  Sigma.mm <- Sigma[missing.indeces,missing.indeces]
  Sigma.oo <- Sigma[observed.indeces,observed.indeces]
  Sigma.mo <- Sigma[missing.indeces,observed.indeces]
  Sigma.om <- t(Sigma.mo)
  
  # The spectral decomposition of Sigma.oo = P D P', eliminating the columns corresponding to null eigenvalues is:
  Sigma.oo.D <- eigen(Sigma.oo,symmetric=TRUE,only.values=TRUE)$values
  number.redundant <- sum((Sigma.oo.D)^2 < 0.00000000001)
  dim.sig <- length(Sigma.oo.D) - number.redundant
  Sigma.oo.D <- Sigma.oo.D[1:dim.sig]
  Sigma.oo.P <- svd(Sigma.oo,nu=dim.sig,nv=dim.sig)$u
  
  # The Moore-Penrose inverse of Sigma.oo is:
  Sigma.oo.M.P.Inverse <- Sigma.oo.P %*% diag(1.0/Sigma.oo.D) %*% t(Sigma.oo.P)
  
  # The covariance matrix of the missing data given the observed data (and mu, sigma2.i, and W.i) is:
  Omega.m <- Sigma.mm - Sigma.mo %*% Sigma.oo.M.P.Inverse %*% Sigma.om
  
  # Important: For this particular example, the rank of Omega.m is equal to 1.
  # That means that even though there are 4 missing observations, there are 3 hard constrains
  # and only one missing value to be estimated.
  
  # The mean vector of the missing data given the observed data (and mu, sigma2.i, and W.i) is:
  gamma.m <- mu[missing.indeces] + Sigma.mo %*% Sigma.oo.M.P.Inverse %*% (z.star[observed.indeces]-mu[observed.indeces])
  
  
  # First, compute the spectral decomposition of Omega.m, eliminating the columns corresponding to null eigenvalues:
  D <- eigen(Omega.m,symmetric=TRUE,only.values=TRUE)$values
  dimension <- length(D)
  number.null <- sum((D)^2 < 0.00000000001)
  nfq <- dimension-number.null      # Number of free quantities
  D <- D[1:nfq]
  P <- svd(Omega.m,nu=nfq,nv=nfq)$u
  
  # Now simulate the missing data:
  y.star <- z.star[1:(4*N)]
  y.star[missing.indeces] <- as.vector(gamma.m) + P %*% as.matrix(sqrt(D) * rnorm(nfq,0,1))
  
  imputed.vector <- y.star[missing.indeces]
  
  list(y.star=y.star,imputed.vector=imputed.vector)
}


## Function to calculate V matrix given sigmas ---------------------------------
my.V <- function(sigma,N) 
{
  temp <- NULL
  
  for (n in 1:N)
  {
    temp <- c(temp,rep(sigma[n],4))
  }
  
  V <- diag(temp)
  return(V)
}


## Function to transform y.star to y, and mu into x.star, ----------------------
my.y <- function(y.star,N,T)
{
  y <- matrix(0,N,4*T)
  
  for (n in 1:N)
  {
    y[n,] <- as.vector(t(y.star[,(4*n-3):(4*n)]))
  }
  return(y)
}

my.x.star <- function(mu,N,T)
{ 
  x.star <- matrix(0,T,4*N)
  
  for (t in 1:T)
  {
    x.star[t,] <- as.vector(t(mu[,(4*t-3):(4*t)]))
  }
  
  return(x.star)
}


## FFBS by M. A. R. Ferreira, 2008. --------------------------------------------
ffbs.MARF <- function(T,y,sigma2,W,a1,R1)
{ 
  # Define the auxiliary vectors:
  theta <- rep(0,T)
  a <- rep(0,T)
  R <- rep(0,T)
  f <- rep(0,T)
  Q <- rep(0,T)
  A <- rep(0,T)
  e <- rep(0,T)
  m <- rep(0,T)
  C <- rep(0,T)
  # Kalman filter:
  a[1] <- a1
  R[1] <- R1
  f[1] <- a[1]
  Q[1] <- R[1]+sigma2
  A[1] <- R[1] / Q[1]
  e[1] <- y[1] - f[1]
  C[1] <- 1.0 / (1.0/R[1] + 1.0/sigma2)          # Equivalent to C[1] <- R[1] - A[1]^2 * Q[1]  and 
  m[1] <- C[1] * (a[1]/R[1] + y[1]/sigma2)       # m[1] <- a[1] + A[1] * e[1], but numerically more stable
  for(t in 2:T)
  {
    a[t] <- m[t-1]
    R[t] <- C[t-1] + W * sigma2
    f[t] <- a[t]
    Q[t] <- R[t]+sigma2
    A[t] <- R[t] / Q[t]
    e[t] <- y[t] - f[t]
    m[t] <- a[t] + A[t] * e[t]
    C[t] <- R[t] - A[t]^2 * Q[t]
  }
  
  # Now, the backward sampler:
  theta[T] <- rnorm(1,m[T],sqrt(C[T]))
  for(t in (T-1):1)
  {
    H <- 1.0 / ( 1.0/(W * sigma2) + 1.0/C[t] )
    h <- H * ( theta[t+1]/(W * sigma2) + m[t]/C[t] )
    theta[t] <- rnorm(1,h,sqrt(H))
  }
  theta
}


## Quantile functions ----------------------------------------------------------
quant025 = function(x){quantile(x,0.025)}
quant05 = function(x){quantile(x,0.05)}
quant95 = function(x){quantile(x,0.95)}
quant975 = function(x){quantile(x,0.975)}


## Function for reading in data ------------------------------------------------
read_data = function(df){
  N = ncol(df$detail)
  T = nrow(df$a_total)
  
  y.plot = (cbind(df$q_total, df$detail) %>%
    as.matrix(., 4*T, N + 1) %>% 
    t)/1000
  
  y.plot.filled <- y.plot[-1,]
  
  for (n in 1:N){
    y.plot.filled[n,is.na(y.plot.filled[n,])] <- mean(y.plot.filled[n,],na.rm=TRUE)
  }
  
  z.star    = matrix(0,T,(5*N+4))
  y.star    = matrix(0,T,4*N)
  y.star.obs= list(NULL)
  index.mis = list(NULL)
  index.obs = list(NULL)
  length(y.star.obs) = T
  length(index.mis)  = T
  length(index.obs)  = T
  
  for(t in 1:T){	  
    z.star[t,]	<- c(as.vector(as.matrix(df$detail[(4*t-3):(4*t),])),
                    df$q_total[(4*t-3):(4*t),],
                    as.matrix(df$a_total[t,]))/1000
    y.star[t,] 	<- as.vector(t(y.plot.filled[,(4*t-3):(4*t)]))
    index.mis[[t]]  <- which(is.na(z.star[t,][1:(4*N)]))
    index.obs[[t]]  <- which(complete.cases(z.star[t,][1:(4*N)]))
    y.star.obs[[t]] <- y.star[t,][index.obs[[t]]] 
  }
  
  output = list(y.plot=y.plot,z.star=z.star,y.star=y.star,
                y.star.obs=y.star.obs,index.mis=index.mis,
                index.obs=index.obs,T=T,N=N)

  return(output)
}

## BMMI ------------------------------------------------------------------------
BMMI = function(data_raw, n_imputation = 10,
                n.iter = 1400, n.burn = 499){
  # parse data
  data = read_data(data_raw)
  
  y.plot <- data[[1]]
  z.star <- data[[2]]
  y.star <- data[[3]]
  y.star.obs <- data[[4]]
  index.mis  <- data[[5]]
  index.obs  <- data[[6]]
  T          <- data[[7]]
  N          <- data[[8]]
  
  # Starting Values and MCMC Parameters
  sigma2 <- rep(1000,N)
  Av <- 0.01  
  Bv <- 0.01
  Aw <- 3
  Bw <- .1
  
  tsi <- rep(0.1,N)
  
  W <- tsi*sigma2
  
  param <- NULL
  temp <- NULL
  imputations = NULL
  store.param <- NULL
  n.quarters <- 4*T
  
  store.param <- matrix(
    NA,
    nrow=n.iter,
    ncol=N*4*T+sum(lengths(index.mis))+2*N
  )
  
  # The matrix that operates in the individual observations and 
  # returns the individual observations and the several totals is
  H <- rbind(Matrix::Diagonal(4*N) %>% as.matrix,
             t(do.call("rbind", rep(list(Matrix::Diagonal(4)), N)) %>% as.matrix),
             kronecker(Matrix::Diagonal(N) %>% as.matrix,t(rep(1,4))))
  HH <- list(NULL)
  length(HH) <- N
  
  for (t in 1:T)
  {
    annual.missing <- which(is.na(data[[2]][t,][(4*(N+1)+1):(4*(N+1)+N)]))
    if (length(annual.missing)>0) HH[[t]] <- H[-(annual.missing+4*(N+1)),]
    else HH[[t]] <- H
  }
  
  indicator <- rep(0,T)
  for (t in 1:T)
  {
    if (length(index.mis[[t]])>0) indicator[t] <- 1
  }
  
  y.star.mis <- list(NULL)
  length(y.star.mis) <- T
  
  impute <- list(NULL)
  length(impute) <- T
  
  y.star.mis.reg <- list(NULL)
  length(y.star.mis.reg) <- T
  
  mu.reg <- list(NULL)
  length(mu.reg) <- N
  
  for (i in 1:n.iter) {
    # Transform y1.star,...,y6.star into 3 series 
    # y1,y2 and y3 these are series for each sector over the 
    # whole time period.  
    y.transf <- NULL
    y.transf <- my.y(y.star,N,T)
    y <- y.transf
    
    # Sample mu
    mu <- matrix(0,N,4*T)
    for (n in 1:N)
    {
      mu[n,] <-  ffbs.MARF(4*T,y[n,],sigma2[n],tsi[n],0,10^10) 
    }
    x.star.transf <- NULL
    x.star.transf <- my.x.star(mu,N,T)
    x.star <- x.star.transf
    
    # Sample sigma2
    alpha.v  <- rep(NA,N)
    delta.v2 <- rep(NA,N)
    beta.v   <- rep(NA,N)
    sigma2   <- rep(NA,N)
    for (n in 1:N)
    {
      alpha.v[n] <- Av + (2*n.quarters-1)/2
      delta.v2[n]<- .5*sum((y[n] - mu[n,])^2) + .5*sum((mu[n,][n.quarters:2]-mu[n,][(n.quarters-1):1])^2/tsi[n])
      beta.v[n]  <- Bv + delta.v2[n]
      sigma2[n]  <- MCMCpack::rinvgamma(1,alpha.v[n],beta.v[n])
    }
    V <- my.V(sigma2,N)
    
    # Sample tsi
    alpha.w      <- rep(NA,N)
    sum.delta.w2 <- rep(NA,N)
    beta.w       <- rep(NA,N)
    tsi          <- rep(NA,N)
    for (n in 1:N)
    {
      alpha.w[n]      <- Aw+((n.quarters-1)/2)
      sum.delta.w2[n] <- sum((mu[n,][n.quarters:2]-mu[n,][(n.quarters-1):1])^2)/sigma2[n]
      beta.w[n]       <- Bw + .5*sum.delta.w2[n]
      tsi[n]          <- MCMCpack::rinvgamma(1,alpha.w[n],beta.w[n])
    }
    W <- tsi*sigma2
    
    # Impute Missing Data Values
    for (t in 1:T)
    {
      if (indicator[t]>0)
      {
        y.star.temp <- NULL
        y.star.mis.temp <- NULL
        impute.temp <- NULL
        y.star.mis.reg.temp <- NULL
        
        impute.temp <- my.missing(z.star[t,],x.star[t,],V,index.mis[[t]],HH[[t]],N)
        
        y.star.temp <- impute.temp$y.star
        y.star.mis.temp <- impute.temp$imputed.vector
        
        y.star.mis.reg.temp <- y.star.mis.temp*1000
        
        y.star[t,] <- y.star.temp
        y.star.mis[[t]] <- y.star.mis.temp
        impute [[t]] <- impute.temp
        y.star.mis.reg[[t]] <- y.star.mis.reg.temp
      }
    }
    store.param[i,] <- c(as.vector(t(mu*1000)), unlist(y.star.mis), tsi, sigma2)
    imputations[[i]] = y.star.mis
  }
  
  # finish
  result = NULL
  interval = (n.iter - n.burn) %/% n_imputation
  for (i in 1:n_imputation) {
    imputation = imputations[[n.burn + (i - 1)*interval]]
    temp = data_raw$detail
    for (j in 1:T) {
      temp[(j*4-3):(j*4),][is.na(temp[(j*4-3):(j*4),])] = imputation[[j]]*1000
    }
    result[[paste0("Imputation-", i)]] = temp
  }
  return(result)
}

