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

##################################
### some data generating processes
##################################

# Ornstein–Uhlenbeck process
error_sim <- function(grid_point, samplesize){
    ti=1:grid_point/grid_point
    comat=matrix(NA,grid_point,grid_point)
    for (i in 1:grid_point){
      comat[i,]=exp(1)^(-ti[i]/2-ti/2)*pmin(exp(1)^(ti[i]),exp(1)^(ti))
    }
    epsilon=mvrnorm(n = samplesize, mu = c(rep(0,grid_point)), Sigma = comat, empirical = FALSE)

    return(t(epsilon))
  }
  
int_approx <- function(x){
        temp_n = NROW(x)
        return((1/(temp_n)) * sum(x))
    }

### functional AR(1)
# @parameter c - constant
# @parameter a - alpha under h0
# @parameter a1 - alpha under h1
# @parameter ks - location of change point
far1_simSB <- function(point_grid, size_sample, c, a0, a1, ks){
    times = 1:point_grid/point_grid
    n0 = 1000+size_sample

    sim_ar_matrix0 = sim_ar_matrix1 = matrix(NA, point_grid, n0)
    error_matrix = error_sim(grid_point = point_grid, samplesize = n0)

    alpha_par0 = function(t,s){
      return(a0 * t * (1-t) * s * (1-s))
    }
    alpha_par1 = function(t,s){
      return((a0+a1) * t * (1-t) * s * (1-s))
    }
    delta_par = function(t){
      return(c* t * (1-t))
    }
    sim_ar_matrix0[,1] = sim_ar_matrix1[,1] = delta_par(times) + error_matrix[,1]
    
    # H0
    for(j in 2:n0){
        for(i in 1:point_grid){
            vector_for_alpha_op = alpha_par0(times[i], times) * sim_ar_matrix0[,(j-1)]
            sim_ar_matrix0[i,j] = delta_par(times[i]) + int_approx(vector_for_alpha_op) + error_matrix[i,j]
        }
    }
    # H1
    for(j in 2:n0){
        for(i in 1:point_grid){
            if (j <= 1000 + floor(size_sample*ks)){
               vector_for_alpha_op = alpha_par0(times[i], times) * sim_ar_matrix1[,(j-1)]
               sim_ar_matrix1[i,j] = delta_par(times[i]) + int_approx(vector_for_alpha_op) + error_matrix[i,j]
            }
            if (j > 1000 + floor(size_sample*ks)){
               vector_for_alpha_op = alpha_par1(times[i], times) * sim_ar_matrix1[,(j-1)]
               sim_ar_matrix1[i,j] = delta_par(times[i]) + int_approx(vector_for_alpha_op) + error_matrix[i,j]
            }
        }
    }

    return(list(sim_ar_matrix0[,1001:(1000 + size_sample)], sim_ar_matrix1[,1001:(1000 + size_sample)]))
}

### linear process
#similar to FAR

### FGARCH(1,1) only covariance operator
## constant change 
fgarch11_simSBc <- function(point_grid, size_sample, c0, c1, ks){
    times = 1:point_grid/point_grid

    n0 = 1000 + size_sample
    sim_sigma2_matrix = sim_garch_matrix0 = sim_garch_matrix1 = matrix(NA, point_grid, n0)
    error_matrix = error_sim(grid_point = point_grid, samplesize = n0)
    
    alpha_par = function(t,s){
      return(6 * t * (1-t) * s * (1-s))
    } 
    beta_par = function(t,s){
      return(14 * t * (1-t) * s * (1-s))
    }

    delta_par0 = function(t){
      return(c0* t * (1-t))
    }
    delta_par1 = function(t){
      return((c0+c1)* t * (1-t))
    }

    sim_sigma2_matrix[,1] = delta_par0(times)
    sim_garch_matrix0[,1] = sim_garch_matrix1[,1] =delta_par0(times) * error_matrix[,1]

    # H0
    for(j in 2:n0){        
        for(i in 1:point_grid){
            vector_for_alpha_op = alpha_par(times[i], times) * ((sim_garch_matrix0[,(j-1)])^2)
            vector_for_beta_op = beta_par(times[i], times) * ((sim_sigma2_matrix[,j-1]))
            sim_sigma2_matrix[i,j] = delta_par0(times[i]) + int_approx(vector_for_alpha_op) + int_approx(vector_for_beta_op)
        }   
        sim_garch_matrix0[,j] = sqrt(sim_sigma2_matrix[,j]) * error_matrix[,j]
    }
    # H1
    for(j in 2:n0){        
        for(i in 1:point_grid){
          if (j <= 1000 + floor(size_sample*ks)){
            vector_for_alpha_op = alpha_par(times[i], times) * ((sim_garch_matrix1[,(j-1)])^2)
            vector_for_beta_op = beta_par(times[i], times) * ((sim_sigma2_matrix[,j-1]))
            sim_sigma2_matrix[i,j] = delta_par0(times[i]) + int_approx(vector_for_alpha_op) + int_approx(vector_for_beta_op)
          }
          if (j > 1000 + floor(size_sample*ks)){
            vector_for_alpha_op = alpha_par(times[i], times) * ((sim_garch_matrix1[,(j-1)])^2)
            vector_for_beta_op = beta_par(times[i], times) * ((sim_sigma2_matrix[,j-1]))
            sim_sigma2_matrix[i,j] = delta_par1(times[i]) + int_approx(vector_for_alpha_op) + int_approx(vector_for_beta_op)
          }  
        }   
        sim_garch_matrix1[,j] = sqrt(sim_sigma2_matrix[,j]) * error_matrix[,j]
    }

    return(list(sim_garch_matrix0[,1001:(1000 + size_sample)], 
                sim_garch_matrix1[,1001:(1000 + size_sample)]))
}


## coefficient change
fgarch11_simSBab <- function(point_grid, size_sample, a0, a1, b0, b1, ks){
    times = 1:point_grid/point_grid

    n0 = 1000 + size_sample
    sim_sigma2_matrix = sim_garch_matrix0 = sim_garch_matrix1 = matrix(NA, point_grid, n0)
    error_matrix = error_sim(grid_point = point_grid, samplesize = n0)
    
    alpha_par0 = function(t,s){
      return(a0 * t * (1-t) * s * (1-s))
    } 
    beta_par0 = function(t,s){
      return(b0 * t * (1-t) * s * (1-s))
    }
    alpha_par1 = function(t,s){
      return((a0+a1) * t * (1-t) * s * (1-s))
    } 
    beta_par1 = function(t,s){
      return((b0+b1) * t * (1-t) * s * (1-s))
    }

    delta_par = function(t){
      return(0.1* t * (1-t))
    }

    sim_sigma2_matrix[,1] = delta_par(times)
    sim_garch_matrix0[,1] = sim_garch_matrix1[,1] =delta_par(times) * error_matrix[,1]

    # H0
    for(j in 2:n0){        
        for(i in 1:point_grid){
            vector_for_alpha_op = alpha_par0(times[i], times) * ((sim_garch_matrix0[,(j-1)])^2)
            vector_for_beta_op = beta_par0(times[i], times) * ((sim_sigma2_matrix[,j-1]))
            sim_sigma2_matrix[i,j] = delta_par(times[i]) + int_approx(vector_for_alpha_op) + int_approx(vector_for_beta_op)
        }   
        sim_garch_matrix0[,j] = sqrt(sim_sigma2_matrix[,j]) * error_matrix[,j]
    }
    # H1
    for(j in 2:n0){        
        for(i in 1:point_grid){
          if (j <= 1000 + floor(size_sample*ks)){
            vector_for_alpha_op = alpha_par0(times[i], times) * ((sim_garch_matrix1[,(j-1)])^2)
            vector_for_beta_op = beta_par0(times[i], times) * ((sim_sigma2_matrix[,j-1]))
            sim_sigma2_matrix[i,j] = delta_par(times[i]) + int_approx(vector_for_alpha_op) + int_approx(vector_for_beta_op)
          }
          if (j > 1000 + floor(size_sample*ks)){
            vector_for_alpha_op = alpha_par1(times[i], times) * ((sim_garch_matrix1[,(j-1)])^2)
            vector_for_beta_op = beta_par1(times[i], times) * ((sim_sigma2_matrix[,j-1]))
            sim_sigma2_matrix[i,j] = delta_par(times[i]) + int_approx(vector_for_alpha_op) + int_approx(vector_for_beta_op)
          }  
        }   
        sim_garch_matrix1[,j] = sqrt(sim_sigma2_matrix[,j]) * error_matrix[,j]
    }

    return(list(sim_garch_matrix0[,1001:(1000 + size_sample)], 
                sim_garch_matrix1[,1001:(1000 + size_sample)]))
}