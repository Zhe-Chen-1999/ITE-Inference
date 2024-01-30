################################################################################
############################# helper functions #################################
################################################################################
library(RIQITE)
library(extraDistr)
library(dplyr)

rank_score <- function( n, method.list = list(name = "Wilcoxon") ){
  
  if ( method.list$name == "DIM" ) {
    stop( "Can't calculate rank scores for DIM" )
  }
  
  if(method.list$name == "Wilcoxon"){
    score = c(1:n)
    score = score/max(score)
    return(score)
  }
  
  if(method.list$name == "Stephenson"){
    score = choose( c(1:n) - 1, method.list$s - 1 )
    score = score/max(score)
    return(score)
  }
  
}


#### (1-alpha) simultaneous confidence/prediction intervals (lower bound) for all effects quantiles among treated, control or all units ####
ci_lower_quantile_exp <- function( Z, Y, treat.method.list = list(name = "Stephenson", s = 6), control.method.list = list(name = "Stephenson", s = 6),
                                   score = NULL, stat.null = NULL, nperm = 10^4, Z.perm = NULL, alpha=0.05,  set = "treat", alpha.ratio.treat = 0.5,  tol = 10^(-3) ){
  n = length(Z)
  
  if(set == "treat"){
    ci.treat = RIQITE::ci_quantile(Z = Z, Y = Y, k.vec = (n-sum(Z)+1):n, alternative = "greater", method.list = treat.method.list, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm, alpha = alpha, tol = tol, switch = FALSE)
    ci.treat = as.numeric( ci.treat$lower )
    return(ci.treat)
  }
  
  if(set == "control"){
    if(is.null(stat.null)){
      stat.null.control = NULL
    }else{
      if(is.null(score)){
        score = rank_score(n, method.list = control.method.list) #### update later
      }
      stat.null.control = sum(score) - stat.null
    }
    
    if(is.null(Z.perm)){
      Z.perm.control = NULL
    }else{
      Z.perm.control = 1 - Z.perm
    }
    ci.control = RIQITE::ci_quantile(Z = 1-Z, Y = -Y, k.vec = (n-sum(1-Z)+1):n, alternative = "greater", method.list = control.method.list, score = score, stat.null = stat.null.control, nperm = nperm, Z.perm = Z.perm.control, alpha = alpha, tol = tol, switch = FALSE)
    ci.control = as.numeric( ci.control$lower )
    return(ci.control)
  }
  
  if(set == "all"){
    
    ci.treat = RIQITE::ci_quantile(Z = Z, Y = Y, k.vec = (n-sum(Z)+1):n, alternative = "greater", method.list = treat.method.list, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm, alpha = alpha*alpha.ratio.treat, tol = tol, switch = FALSE)
    ci.treat = as.numeric( ci.treat$lower )
    
    
    # cat("treated","\n")
    # print(ci.treat)
    
    if(is.null(stat.null)){
      stat.null.control = NULL
    }else{
      if(is.null(score)){
        score = rank_score(n, method.list = control.method.list) #### update later
      }
      stat.null.control = sum(score) - stat.null
    }
    
    if(is.null(Z.perm)){
      Z.perm.control = NULL
    }else{
      Z.perm.control = 1 - Z.perm
    }
    ci.control = RIQITE::ci_quantile(Z = 1-Z, Y = -Y, k.vec = (n-sum(1-Z)+1):n, alternative = "greater", method.list = control.method.list, score = score, stat.null = stat.null.control, nperm = nperm, Z.perm = Z.perm.control, alpha = alpha*(1-alpha.ratio.treat), tol = tol, switch = FALSE)
    ci.control = as.numeric( ci.control$lower )
    
    ci.all = sort( c(ci.treat, ci.control) )
    
    # cat("control","\n")
    # print(ci.control)
    
    return( ci.all)
  }
}


#### correction term based on multivariate hypergeometric distribution ####

correct_hypergeom <- function( N = 100, n = 50, K.vec = c(60, 80), kk.vec = c(30, 40), ndraw = 10^4, hg.draw = NULL ){
  if( length(K.vec) > 1 ){
    if(is.null(hg.draw)){
      hg.draw = rmvhyper(ndraw, diff( c(0, K.vec, N) ), n)
      hg.draw = hg.draw[, -1, drop=FALSE]
    }
    H.sum = apply( hg.draw[, ncol(hg.draw):1, drop=FALSE], 1, cumsum )
    H.sum = matrix(H.sum, nrow = ncol(hg.draw) )
    H.sum = H.sum[nrow(H.sum):1, , drop=FALSE]
    prob = mean( apply( H.sum - (n - kk.vec) > 0, 2, max ) )
    return(prob)
  }
  
  if(length(K.vec) == 1){
    prob = 1 - phyper(n-kk.vec, N-K.vec, K.vec, n)
    return(prob)
  }
}

## calculate the correction term Delta using the proposed choice of kk.vec (k_prime)
threshold_correct_hypergeom <- function( N = 100, n = 50, K.vec = c(60, 80), alpha = 0.05, ndraw = 10^4, hg.draw = NULL, tol = 10^(-2) ){
  
  if(length(K.vec) > 1){
    if(is.null(hg.draw)){
      hg.draw = rmvhyper(ndraw, diff( c(0, K.vec, N) ), n)
      hg.draw = hg.draw[, -1, drop=FALSE] 
    }
    
    kappa = NULL
    
    kappa.lower = 1/length(K.vec)
    kappa.upper = 1
    prob.lower = correct_hypergeom( N = N, n = n, K.vec = K.vec, kk.vec = n-qhyper(1-kappa.lower*alpha, N - K.vec, K.vec, n), hg.draw = hg.draw )
    prob.upper = correct_hypergeom( N = N, n = n, K.vec = K.vec, kk.vec = n-qhyper(1-kappa.upper*alpha, N - K.vec, K.vec, n), hg.draw = hg.draw )
    
    if(prob.upper <= alpha){
      kappa = kappa.upper
    }
    
    if(prob.lower > alpha){
      kappa = kappa.lower
    }
    
    while(is.null(kappa)){
      kappa.mid = (kappa.upper+kappa.lower)/2
      prob.mid = correct_hypergeom( N = N, n = n, K.vec = K.vec, kk.vec = n-qhyper(1-kappa.mid*alpha, N - K.vec, K.vec, n), hg.draw = hg.draw )
      
      if(prob.mid <= alpha){
        kappa.lower = kappa.mid
      }else{
        kappa.upper = kappa.mid
      }
      
      if( kappa.upper - kappa.lower <= tol ){
        kappa = kappa.lower 
      }
    }
    
    kk.vec = n-qhyper(1-kappa*alpha, N - K.vec, K.vec, n)
    prob = correct_hypergeom( N = N, n = n, K.vec = K.vec, kk.vec = kk.vec, hg.draw = hg.draw )
    return(list(kappa = kappa, prob = prob, kk.vec = kk.vec))
  }
  
  if(length(K.vec) == 1){
    kappa = 1
    temp = qhyper(1-alpha, N - K.vec, K.vec, n)
    # cat('this is temp:', temp, 'n')
    #cat('this is n:', n, '\n')
    kk.vec = n-qhyper(1-alpha, N - K.vec, K.vec, n)
    #cat('this is kk.vec:', kk.vec, '\n')
    prob = correct_hypergeom( N = N, n = n, K.vec = K.vec, kk.vec = kk.vec)
    return(list(kappa = kappa, prob = prob, kk.vec = kk.vec))
  }
}


#### generalize ####

ci_lower_quantile_gen <- function( Z, Y, N = 2*length(Z), K.vec = ceiling(c(0.6, 0.7, 0.8, 0.9)*N), 
                                   alpha=0.05, gamma = 0.5, ndraw = 10^4, treat.method.list = list(name = "Stephenson", s = 6),
                                   control.method.list = list(name = "Stephenson", s = 6), score = NULL, stat.null = NULL, nperm = 10^4, Z.perm = NULL,  set = "all", alpha.ratio.treat = 0.5,  tol = 10^(-3) ){
  
  if(set == "all"){
    n.star = length(Z)
  }
  
  if(set == "control"){
    n.star = sum(1-Z)
  }
  
  if(set == "treat"){
    n.star = sum(Z)
  }
  
  if( N > n.star ){
    step1 = threshold_correct_hypergeom(N = N, n = n.star, K.vec = K.vec, alpha = gamma * alpha, ndraw = ndraw)
  }
  
  if(N == n.star){
    step1 = list(prob = 0, kk.vec = K.vec)
  }
  
  step2 = ci_lower_quantile_exp( Z, Y, treat.method.list = treat.method.list, control.method.list = control.method.list, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm, 
                                 alpha= alpha - step1$prob,  set = set, alpha.ratio.treat = alpha.ratio.treat,  tol = tol )
  step2 = c(-Inf, step2) # add lower confidence limit "-Inf" for case kk.vec = 0
 
  ci = step2[step1$kk.vec+1] # index +1 for case kk.vec = 0
  
  K.vec = sort(K.vec)
  conf.int = data.frame(k = K.vec, lower = ci, upper = Inf)
  return(conf.int)
}


# specify individual/simultaneous inference
ci_lower_quantile_generalize <- function(Z, Y, N, k_vec = ceiling(c(0.6, 0.7, 0.8, 0.9)*N), 
                                         alpha=0.05, gamma = 0.5, ndraw = 10^4, 
                                         treat.method.list = list(name = "Stephenson", s = 6), control.method.list = list(name = "Stephenson", s = 6), score = NULL, stat.null = NULL, nperm = 10^4, Z.perm = NULL,  
                                         set = "all", alpha.ratio.treat = 0.5, tol = 10^(-3), simul = TRUE){
  
  k_vec = sort(k_vec)
  if (simul) {
    res = ci_lower_quantile_gen( Z = Z, Y = Y, N = N, K.vec = k_vec, 
                                 alpha = alpha, gamma = gamma, ndraw = ndraw, 
                                 treat.method.list = treat.method.list, control.method.list = control.method.list, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm,
                                 set = set, alpha.ratio.treat = alpha.ratio.treat,  tol = tol )
  }
  else {
    res = NULL
    for (k in k_vec){
        res = rbind(res, 
                    ci_lower_quantile_gen( Z = Z, Y = Y, N = N, K.vec = c(k), 
                                           alpha = alpha, gamma = gamma, ndraw = ndraw, 
                                           treat.method.list = treat.method.list, control.method.list = control.method.list, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm,
                                           set = set, alpha.ratio.treat = alpha.ratio.treat,  tol = tol))
    }
  }
  return(res)
}


# #### example ####
# stat.null = null_dist(n=100, m=n/2, method.list = list(name = "Stephenson", s = 6), nperm = 10^6)
# stat.null = NULL
# 
# n = 100
# Z = sample(c(rep(1,n/2), rep(0,n/2)))
# Y0 = rnorm(n)
# sigma = 5
# Y1 = rnorm(n)*sigma + 2 * sqrt(1+sigma^2)
# Y = Z*Y1 + (1-Z)*Y0
# 
# # M1
# ci.comb.gen = ci_lower_quantile_gen( Z, Y, N = length(Z), K.vec = ceiling(c(0.6, 0.7)*N), 
#                                      alpha=0.05, gamma = 0.5, ndraw = 10^5, 
#                                      method.list = list(name = "Stephenson", s = 6), stat.null = stat.null, nperm = 10^4, 
#                                      set = "all", 
#                                      alpha.ratio.treat = 0.5,  tol = 10^(-3) )
# ci.comb.gen
# 
# # M2
# ci.gen1 = ci_lower_quantile_gen( Z, Y, N = length(Z), K.vec = ceiling(c(0.6, 0.7, 0.8, 0.9)*N), 
#                                  alpha=0.05/2, gamma = 0.5, ndraw = 10^5, 
#                                  method.list = list(name = "Stephenson", s = 6), stat.null = stat.null, nperm = 10^4, 
#                                  set = "treat", 
#                                  alpha.ratio.treat = 0.5,  tol = 10^(-3) )
# 
# ci.gen2 = ci_lower_quantile_gen( Z, Y, N = length(Z), K.vec = ceiling(c(0.6, 0.7, 0.8, 0.9)*N), 
#                                  alpha=0.05/2, gamma = 0.5, ndraw = 10^5, 
#                                  method.list = list(name = "Stephenson", s = 6), stat.null = stat.null, nperm = 10^4, 
#                                  set = "control", 
#                                  alpha.ratio.treat = 0.5,  tol = 10^(-3) )
# 
# pmax( ci.gen1, ci.gen2 )
# 
# # compare
# ci.comb.gen - pmax( ci.gen1, ci.gen2 )

################################ Simulation ####################################
################################################################################
###################### Data generating process ################################# 
################################################################################
# Y0 ~ N(0, sqrt(rho^2)), Y1 ~ N(0, sqrt(1-rho^2))
# The total variance of Y1 and Y0 is 1
# So that the results for different rho^2 may be more comparable

# simulate_dt_mixture <- function(N, n, k, rho, mu = 2){
#   
#   # Generate Y(0) and Y(1) for whole population
#   Y_0 = rnorm(N, mean = 0, sd = rho)
#   Y_1 = rnorm(N, mean = mu, sd = sqrt(1-rho^2))
#   Tau = Y_1 - Y_0
#   
#   # Sample experimental units
#   sample_index = sample(1:N, n)
#   y_0 = Y_0[sample_index]
#   y_1 = Y_1[sample_index]
#   tau = y_1 - y_0
#   
#   # n_t treated units and n_c control units
#   n_t = n * k
#   n_c = n - n_t
#   
#   Z = c(rep(1, n_t), rep(0, n_c))
#   Y_obs = Z * y_1 + (1-Z)*y_0
#   
#   dt_N = data.frame(Y_0 = Y_0, Y_1 = Y_1, Tau = Tau)
#   dt_n = data.frame(Z = Z, y_0 = y_0, y_1 = y_1, Y_obs = Y_obs, tau = tau)
#   dt_n = slice(dt_n, sample(1:n()))
#   
#   return(list(dt_N, dt_n))
# }

## example

# dt = simulate_dt_mixture(100, 50, 0.5, 0.2)
# dt[1]: population data Nx3
# dt[2]: sample data nx5

################################################################################
######################  Methods to be compared ################################# 
################################################################################
# M0: Method in Caughey et al. (2023) 
# M1: Method combining treated and control using Bonferroni correction (apply M0 twice)
# M2: Method using Berger and Boos (1994)'s approach

## The standard output for all methods: 
## J (number of testing quantiles) rows and three columns. Three columns are k, lower, and upper
################################################################################
# M0: Method in Caughey et al. (2023) 
method_caughey <- function(Z, Y, k_vec, method.list = list(name = "Wilcoxon"), nperm = 10^4,  alpha = 0.05){
  return(RIQITE::ci_quantile(Z = Z, Y = Y, k.vec = k_vec,
                             alternative = "greater", 
                             method.list = method.list,
                             nperm = nperm,  alpha = alpha))
}
################################################################################
# M1: Method combining treated and control using Bonferroni correction (apply M0 twice)
method_combine <- function(Z, Y, N, k_vec, control.method.list = list(name = "Stephenson", s = 6), treat.method.list = list(name = "Stephenson", s = 6),
                           simul = TRUE,
                           stat.null = NULL, score = NULL, Z.perm = NULL, 
                           alpha=0.05, gamma = 0.5, alpha.ratio.treat = 0.5, 
                           ndraw = 10^5, nperm = 10^4, tol = 10^(-3)){
  
  return(ci_lower_quantile_generalize(Z = Z, Y = Y, N = N, k_vec = k_vec, simul = simul,
                                      alpha = alpha, gamma = gamma, ndraw = ndraw, control.method.list = control.method.list,
                                      treat.method.list = treat.method.list, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm,
                                      set = "all", alpha.ratio.treat = alpha.ratio.treat, tol = tol))
}

## example

# method_combine( Z, Y, N = length(Z), k_vec = ceiling(c(0.6, 0.7, 0.8, 0.9)*N),simul = TRUE) 
#
################################################################################
# M2: Method using Berger and Boos (1994)'s approach
method_bergerboos <- function(Z, Y, N, k_vec, treat.method.list = list(name = "Stephenson", s = 6), control.method.list = list(name = "Stephenson", s = 6),
                              simul = TRUE,
                              stat.null = NULL, score = NULL, Z.perm = NULL, 
                              alpha=0.05, gamma = 0.5, alpha.ratio.treat = 0.5, 
                              ndraw = 10^5, nperm = 10^4, tol = 10^(-3)){
  
  ci.gen1 = ci_lower_quantile_generalize(Z, Y, N, k_vec, treat.method.list = treat.method.list, 
                                         stat.null = stat.null, set = "treat", simul = simul,
                                         alpha=alpha/2, gamma = gamma, alpha.ratio.treat = alpha.ratio.treat, 
                                         score = score, Z.perm = Z.perm,
                                         ndraw = ndraw, nperm = nperm, tol = tol)
  
   # cat("treat", '\n')
   # print(ci.gen1)
  ci.gen2 = ci_lower_quantile_generalize(Z, Y, N, k_vec, control.method.list = control.method.list, 
                                         stat.null = stat.null, set = "control", simul = simul,
                                         alpha=alpha/2, gamma = gamma, alpha.ratio.treat = alpha.ratio.treat, 
                                         score = score, Z.perm = Z.perm,
                                         ndraw = ndraw, nperm = nperm, tol = tol)
   # cat("control", '\n')
   # print(ci.gen2)
  return(pmax( ci.gen1, ci.gen2 ))
}

## example

# method_bergerboos( Z, Y, N = length(Z), k_vec = ceiling(c(0.6, 0.7, 0.8, 0.9)*N), simul = FALSE,
#                   method.list = list(name = "Stephenson", s = 6), stat.null = stat.null)
#################################################################################################################################

################################## Placebo-controlled trials with highly specific endpoints ################################## 
################################## Method in Sedransk and Meyer (1978, JRSS-B) ##############################################
method_SM <- function(Z, Y, LOD, k_vec, simul = FALSE, nperm = 10^6, Z.perm = NULL, alpha = 0.05, tol = 10^(-3)){
  
  if (simul) {
    res = method_SM_SimulCI(Z = Z, Y = Y, LOD = LOD, k.vec = k_vec, nperm = nperm, Z.perm = Z.perm, alpha = alpha, tol = tol)
  }
  else {
    res = method_SM_singleCI(Z = Z, Y = Y, LOD = LOD, k_vec = k_vec, alpha = alpha)
  }
  return(res)
}

method_SM_SimulCI <- function(Z, Y, LOD, k.vec = NULL, nperm = 10^6, Z.perm = NULL, alpha = 0.05, tol = 10^(-3)){
  Y = Y - LOD
  N = length(Z)
  m = sum(Z)
  
  if(is.null(k.vec)) {
    k.vec = c(1:N)
  }
  
  k.vec = sort(k.vec, decreasing = FALSE)
  J = length(k.vec)
  conf.int = data.frame(k = k.vec, lower = rep(NA, J), upper = rep(NA, J))
  
  # f1 is decreasing with corrected_alpha
  f1 <- function(corrected_alpha){
    return(sci_prob_lb(Z = Z, k.vec = k.vec, nperm = nperm, Z.perm = Z.perm, alpha = corrected_alpha) - (1-alpha))
  }
  
  # level correction of each individual CI
  corrected_alpha = uniroot(f1, interval = c(0, alpha), extendInt = "downX", tol = tol)$root
  
  conf.int = method_SM_singleCI(Z, Y+LOD, LOD, k_vec = k.vec, alpha = corrected_alpha)
  
  return(conf.int)
}

sci_prob_lb <- function(Z, k.vec, nperm = 10^6, Z.perm = NULL, alpha = 0.05){
  N = length(Z)
  m = sum(Z)
  J = length(k.vec)
  
  # get Z.perm if Z.perm is null
  if(is.null(Z.perm)){
    Z.perm = assign_CRE(N, m, nperm)
  } else {
    stopifnot(is.matrix(Z.perm))
  }
  
  stat = matrix(NA, nrow = J,  ncol = ncol(Z.perm))
  G = matrix(NA, nrow = J,  ncol = ncol(Z.perm))
  q = c()
  for(j in 1:J){
    # stat[j, ]: permutated \sum_i{Z_i*I(i>k_j)}
    stat[j, ] = colSums(Z.perm * (c(1:N) > k.vec[j]))
    q[j] = qhyper(1-alpha, m = N - k.vec[j], n = k.vec[j], k = m)
    G[j, ] = stat[j, ] -  q[j]
  }
  
  LB = 1 - mean(colSums(G > 0) > 0)
  return(LB)
}

method_SM_singleCI <- function(Z, Y, LOD, k_vec, alpha = 0.05){
  
  nt = sum(Z)
  nc = sum(1-Z)
  N = nt+nc
  tau = sort(Y[which(Z==1)] - LOD[which(Z==1)])
  all_ci = c()
  for(i in k_vec){
    k = nt - qhyper(1-alpha, m = N - i, n = i, k = nt)
    if(k==0){
      all_ci = c(all_ci, -Inf)
    }else{
      all_ci = c(all_ci, tau[k])
    }
  }
  conf.int = data.frame(k = k_vec, 
                        lower = all_ci, upper = Inf)
  return(conf.int)
  
}
