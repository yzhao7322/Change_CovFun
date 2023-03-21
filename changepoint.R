#packages
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

##### The results in Section 4. Theorem 4.1
################## some useful functions for the limit distribution of estimator
mkfun <- function(kappa,theta,tim){
  if (tim < 0){
    mk = (1-kappa)*(1-theta)+kappa*theta
  }
  else if (tim > 0){
    mk = (1-kappa)*theta+kappa*(1-theta)
  }
  else {
    mk = 0
  }
  return(mk)
}


twowiener <- function(lim_N){
  times1 = seq(-1,0,length=lim_N/2)
  times2 = seq(0,1,length=lim_N/2)
  times = c(times1[1:(lim_N/2-1)], times2)

  w1 = BM(x=0, t0=0, T=1, N=lim_N/2-1)
  w2 = BM(x=0, t0=0, T=1, N=lim_N/2-1)
  w = c(rev(w1[2:(lim_N/2)]), w2)
  return(list(times,w))
}


lim_cp <- function(reps,kappa,theta){
  N = 10000
  cpseq = c(rep(NA,reps))
  for (j in 1:reps){
    wtime = twowiener(N)
    times = wtime[[1]]
    w = wtime[[2]]
    wseq = c(rep(NA,N-1))
    for (i in 1:(N-1)){
      wseq[i] = w[i]-abs(times[i])*mkfun(kappa,theta,times[i])
    }
    mcp = max(wseq)
    limcp = which(wseq == mcp)
    cpseq[j] = limcp
  }
  return(cpseq)
}




###########  estimator


##### estimate size of change
sizechange <- function(xd, kstar){
  N=ncol(xd)

  sample_cov<-function(data){
    N=ncol(data)
    varmtx = 0
    for (i in 1:N){
      varmtx = varmtx + data[,i]%o%data[,i]
    }
    return(varmtx/N)
  }

  error = apply(xd,2,function(x,xmean){x-xmean}, xmean = rowMeans(xd))
  error_before = error[,1:kstar]
  error_after = error[,(kstar+1):N]
  
  var_before = sample_cov(error_before)
  var_after = sample_cov(error_after)
  var_change = var_before - var_after
  return(var_change)
}

########

l2norm <-function(vec){
     return(sqrt(sum(vec^2)))
  }
  
tau_est <-function(xd, kstar, len){
  grid_point = nrow(xd)
  N = ncol(xd)
   
  rref = runif(len, 0, 1)
  rref = c(sort(rref), 1)
  rrefind = round(rref * grid_point)
  rrefind[which(rrefind==0)] = 1
  xdmc = xd[rrefind,]
  
  sample_cov<-function(data){
    N=ncol(data)
    varmtx = 0
    for (i in 1:N){
      varmtx = varmtx + data[,i]%o%data[,i]
    }
    return(varmtx/N)
  }

  error = apply(xdmc,2,function(x,xmean){x-xmean}, xmean = rowMeans(xdmc))
  error_before = error[,1:kstar]
  error_after = error[,(kstar+1):N]
  
  var_before = sample_cov(error_before)
  var_after = sample_cov(error_after)
  var_change = var_before - var_after

  ## change star
  var_1 = var_2 = 0
  
  for (i in 1:kstar){
    var_1 = var_1 + (xdmc[,i]-rowMeans(xdmc))%o%(xdmc[,i]-rowMeans(xdmc))
  }
  var_1 = 1/kstar * var_1

  for (i in (kstar+1):N){
    var_2 = var_2 + (xdmc[,i]-rowMeans(xdmc))%o%(xdmc[,i]-rowMeans(xdmc))
  }
  var_2 = 1/(N-kstar) * var_2

  var_star = (var_1-var_2)/l2norm(var_1-var_2)

  ## longrun cov

  zi = zm = array(0, c((len+1), (len+1), N))
  for (i in 1:N){
     zi[,,i] = error[,i]%o%error[,i]
  }
  
  v_dat=array(0, c(len+1, len+1, N))
  for (i in 1:N){
    if (i <= kstar){
      v_dat[,,i] = zi[,,i]-var_1
    }else{
      v_dat[,,i] = zi[,,i]-var_2
    }
  }

  int_approx_tensor<-function(x){# x is a 4-dimensional tensor
      dt=length(dim(x))
      temp_n=nrow(x)
      return((1/temp_n)^dt * sum(x))}

  longd = long_run_covariance_4tensor(v_dat)

  frontvs = rearvs = 0
  for (i in 1:21){
    for (j in 1:21){
      frontvs = frontvs + var_star%o%longd[i,,j,]
    }
  }
  for (i in 1:21){
    for (j in 1:21){
      rearvs = rearvs + frontvs[,i,,j]%o%var_star
    }
  }
  tau = int_approx_tensor(rearvs)

  return(list(var_change,tau))
}



######### a simple experiment
# need functions from dgp.R and stats.R

set.seed(3)

kappa = 1/4
lens = 20
truek = 0.9
reps = 200
stat_vec = c(rep(NA,reps))
samplesize=500

for (i in 1:reps){
   dgp = fgarch11_simSBc(50, samplesize, 0.1, 0.3, truek)
   dh1 = dgp[[2]]
   stat_d0 = weight_TNstat(dh1,kappa)
   cv_d0 = weight_criticalvalueMC(dh1,len = lens,kappa)

   if (stat_d0[[1]]> cv_d0[2]){
     kstar = stat_d0[[2]]
   
     changetau = tau_est(dh1, kstar, len=20)
     cbar = changetau[[1]]
     tau = changetau[[2]]
     cpstat=l2norm(cbar)^2/tau*((kstar/samplesize)-truek)
     stat_vec[i]=cpstat
   }else{
     stat_vec[i]=NA
   }
}