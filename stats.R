# packages
require(data.table)
require(MASS)
require(fda)
require(tictoc)
require(rainbow)
require(sde)
require(xtable)
require(mvtnorm)
require(tseries)
require(expm)
require(tensorA)

### Theorem 2.1/2.2
# xdm - finite realization of DEMEANed functional time series data, where curves are stored in columns.
# u - a fraction index over the interval [0, 1]
ZNstat <- function(xdm, u){
  grid_point = nrow(xdm)
  N = ncol(xdm)
  k = floor(N*u)
  prek = matrix(rowSums(apply(as.matrix(xdm[,1:k]),2,function(x){x%o%x})),grid_point,grid_point)
  fullk = matrix(rowSums(apply(as.matrix(xdm),2,function(x){x%o%x})),grid_point,grid_point)
  ZNu = N^(-1/2) * (prek - (k/N)*fullk)
  return(ZNu)
}

ZNstat_cp <- function(xdm, u){
  grid_point = nrow(xdm)
  N = ncol(xdm)
  k = floor(N*u)
  prek = matrix(rowSums(apply(as.matrix(xdm[,1:k]),2,function(x){x%o%x})),grid_point,grid_point)
  fullk = matrix(rowSums(apply(as.matrix(xdm),2,function(x){x%o%x})),grid_point,grid_point)
  ZNu = (prek - (k/N)*fullk)
  return(ZNu)
}


# T_N statistic introduced after Theorem 2.1
# xf - finite realization of functional time series data, where curves are stored in columns.
TNstat <- function(xf){
  int_approx=function(x){
      temp_n=NROW(x)
      return((1/temp_n)*sum(x))}	
  grid_point = nrow(xf)
  N = ncol(xf)
  xdm = apply(xf,2,function(x,xmean){x-xmean}, xmean = rowMeans(xf))
  uind = seq(0,1,length = N+1)[2:(N+1)]; zn2=list()
  for (i in 1:N){
  	zn2[[i]] = (ZNstat(xdm, uind[i]))^2
  }
  inm = Reduce(`+`, zn2)/N
  return((1/grid_point)^2*sum(inm))
}


# T_N(\kappa) statistic introduced after Theorem 2.3
weight_TNstat <- function(xf,kappa){
  int_approx_tensor<-function(x){# x is a 4-dimensional tensor
      dt=length(dim(x))
      temp_n=nrow(x)
      return(sum(x)/(temp_n^dt))}

  grid_point = nrow(xf)
  N = ncol(xf)
  xdm = apply(xf,2,function(x,xmean){x-xmean}, xmean = rowMeans(xf))
  uind = seq(0,1,length = N+1)[2:(N+1)]; 
  zn2 = list(); zn_cp = c(rep(0,N))
  for (i in 1:(N-1)){
    zn2[[i]] = (ZNstat(xdm, uind[i]))^2 / ((uind[i]*(1-uind[i]))^(2*kappa))### kappa = 1/4
  
    zn_cp[i] = (N/(i*(N-i)))^(kappa) *  int_approx_tensor( (ZNstat_cp(xdm, uind[i]))^2 )
  }
  inm = Reduce(`+`, zn2)/N
  stat = (1/grid_point)^2*sum(inm)

  mcp = max(zn_cp[ (0.1*N):(0.9*N)])
  changepoint = which(zn_cp == mcp)

  return(list (stat, changepoint) )
}


## Critical values

## useful functions for computing critical values
## long-run covariance operator
# dat - an array with dimension (grid_point,grid_point,N)
long_run_covariance_4tensor <- function (dat)
  {
      grid_point = dim(dat)[1]
      T = dim(dat)[3]
      datmean = apply(dat, c(1,2), mean)
      center_dat = sweep(dat, 1:2, datmean)
      
      cov_l <- function(band, nval) {
          cov_sum = gamma_l(0, nval)

          for (ik in 1:(nval - 1)) {
              cov_sum = cov_sum + kweights(ik/band, kernel = "Bartlett") * (2*gamma_l(ik, nval))# + gamma_l(ik,nval))    ##aperm(gamma_l(ik,nval),c(2,1,3,4)))
          }
          return(cov_sum)
      }
      
      gamma_l <- function(lag, T) {
          gamma_lag_sum = 0
          if (lag >= 0) {
              for (ij in 1:(T - lag)) {
                  gamma_lag_sum = gamma_lag_sum + center_dat[,,ij] %o% center_dat[,,(ij + lag)]
              }
          }
          else {
              for (ij in 1:(T + lag)) {
                  gamma_lag_sum = gamma_lag_sum + center_dat[,,(ij - lag)] %o% center_dat[, ij]
              }
          }
          return(gamma_lag_sum/(T-lag))
      }
      hat_h_opt = T^(1/4)
      lr_covop = cov_l(band = hat_h_opt, nval = T)

      return(lr_covop)
  }

kweights <- function (x, kernel = c("Truncated", "Bartlett", "Parzen", "Tukey-Hanning", 
    "Quadratic Spectral"), normalize = FALSE) 
{
    kernel <- match.arg(kernel)
    if (normalize) {
        ca <- switch(kernel, Truncated = 2, Bartlett = 2/3, Parzen = 0.539285, 
            `Tukey-Hanning` = 3/4, `Quadratic Spectral` = 1)
    }
    else ca <- 1
    switch(kernel, Truncated = {
        ifelse(ca * x > 1, 0, 1)
    }, Bartlett = {
        ifelse(ca * x > 1, 0, 1 - abs(ca * x))
    }, Parzen = {
        ifelse(ca * x > 1, 0, ifelse(ca * x < 0.5, 1 - 6 * (ca * 
            x)^2 + 6 * abs(ca * x)^3, 2 * (1 - abs(ca * x))^3))
    }, `Tukey-Hanning` = {
        ifelse(ca * x > 1, 0, (1 + cos(pi * ca * x))/2)
    }, `Quadratic Spectral` = {
        y <- 6 * pi * x/5
        ifelse(x < 1e-04, 1, 3 * (1/y)^2 * (sin(y)/y - cos(y)))
    })
}  

# compute critical values for TNstat (T_N)
criticalvalueMC <-function(xf,len){
   grid_point = nrow(xf)
   N = ncol(xf)
   
   rref = runif(len, 0, 1)
   rref = c(sort(rref), 1)
   rrefind = round(rref * dim(xf)[1])
   rrefind[which(rrefind==0)] = 1
   xfMC = xf[rrefind,]

   xdm = apply(xfMC,2,function(x,xmean){x-xmean}, xmean = rowMeans(xfMC))
   zi = zm = array(0, c((len+1), (len+1), N))
   for (i in 1:N){
     zi[,,i] = xdm[,i]%o%xdm[,i]
   }
   zimean = apply(zi,c(1,2),mean)
   for (i in 1:N){
     zm[,,i] = zi[,,i]-zimean
   }
   lrcov = long_run_covariance_4tensor(zm) ##23.883 sec elapsed
   lrcov = as.tensor(lrcov/(len+1)^2)
   eigvals=svd.tensor(lrcov,c(3,4),by="e")
   eigmat=as.vector(eigvals$d)

  lim_sum = 0
   for (ell in 1:length(eigmat)){
    klim = 0
    for (k in 1:1000){
      Nm=rnorm(2000,mean=0,sd=1)
      klim = klim + eigmat[ell]/((pi*k)^2)*Nm^2
    }
    lim_sum=lim_sum+klim
  }

  #lim_sum= rowSums(apply(matrix(seq(1,length(eigmat),1),1),2, function(x){ frac = eigmat[x]/((pi*seq(1,k,1))^2);
  #  rowSums(apply(matrix(seq(1,k,1),1),2,function(xx){frac[xx]*rnorm(5000,mean=0,sd=1)^2}))} ) )
   #klim = rowSums(t(frac*t(munor)))
  cv=quantile(lim_sum,probs=c(0.90,0.95,0.99))  
  return(cv)
}

# compute critical values for weight_TNstat ( T_N(\kappa) )

weight_criticalvalueMC <-function(xf,len,kappa){
   grid_point = nrow(xf)
   N = ncol(xf)

   ## cov weight function
   times = 1:grid_point/grid_point
   wmat = matrix(NA,grid_point-2,grid_point-2)
   for (i in 2:(grid_point-1)){
     for (j in 2:(grid_point-1)){
       wmat[i-1,j-1]= (min(times[i],times[j])-times[i]*times[j])/( (times[i]*(1-times[i]))^kappa * (times[j]*(1-times[j]))^kappa )
     }
   }
   weig = as.vector(svd(wmat/grid_point)$d)

   ## cov operators
   rref = runif(len, 0, 1)
   rref = c(sort(rref), 1)
   rrefind = round(rref * dim(xf)[1])
   rrefind[which(rrefind==0)] = 1
   xfMC = xf[rrefind,]

   xdm = apply(xfMC,2,function(x,xmean){x-xmean}, xmean = rowMeans(xfMC))
   zi = zm = array(0, c((len+1), (len+1), N))
   for (i in 1:N){
     zi[,,i] = xdm[,i]%o%xdm[,i]
   }
   zimean = apply(zi,c(1,2),mean)
   for (i in 1:N){
     zm[,,i] = zi[,,i]-zimean
   }

   lrcov = long_run_covariance_4tensor(zm)
   lrcov = as.tensor(lrcov/(len+1)^2)
   eigvals=svd.tensor(lrcov,c(3,4),by="e")
   eigmat=as.vector(eigvals$d)

   lim_sum = 0
   for (ell in 1:length(eigmat)){
    klim = 0
    for (k in 1:length(weig)){
      Nm=rnorm(2000,mean=0,sd=1)
      klim = klim + eigmat[ell]*weig[k]*Nm^2
    }
    lim_sum=lim_sum+klim
  }

  cv=quantile(lim_sum,probs=c(0.90,0.95,0.99))  
  return(cv)
}

