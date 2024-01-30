source("~/helper_functions.R")

library(RIQITE)
library(dplyr)

################# Methods to be compared ####################
# M1: Caughey et al. (2023)
# M2: Combine intervals for treated and control group 
# M3: Berger and Boos's approach

simulate_from_realdata <- function(realdata, antigen_of_interest = 'Con 6 gp120/B', arm.treated, arm.control, nt, nc){
  
  data_analysis = realdata %>%
    filter(antigen == antigen_of_interest) %>%
    filter(rx_code == arm.treated | rx_code == arm.control) 
  
  n = nt + nc
  Y_1 = sample(data_analysis$outcome[which(data_analysis$rx_code == arm.treated)], n, replace = TRUE)+ rnorm(n, mean = 0, sd = 0.15)
  Y_0 = sample(data_analysis$outcome[which(data_analysis$rx_code == arm.control)], n, replace = TRUE)+ rnorm(n, mean = 0, sd = 0.15)
  
  Z = c(rep(1, nt), rep(0, nc))
  
  dt = data.frame(Z = Z, Y_0 = Y_0, Y_1 = Y_1, Y_obs = Z * Y_1 + (1-Z)*Y_0, tau = Y_1 - Y_0)
  dt = slice(dt, sample(1:n()))
  
  return(dt)
}

compare <- function(dt, k_vec, simul = FALSE){
  
  # n is sample size
  n = dim(dt)[1]
  
  # J is length of k_vec
  J = length(k_vec)
  
  # Method M1 is the same for simul = TRUE or FALSE
  m1_s2 = method_caughey(dt$Z, dt$Y_obs, k_vec = k_vec, method.list = list(name = "Stephenson", s = 2))
  m1_s6 = method_caughey(dt$Z, dt$Y_obs, k_vec = k_vec, method.list = list(name = "Stephenson", s = 6))
  
  # Method M2 is the same for simul = TRUE or FALSE
  m2_treat.s2_control.s6 = method_combine(dt$Z, dt$Y_obs, N = n, k_vec = k_vec, simul = simul,
                                           treat.method.list = list(name = "Stephenson", s = 2), 
                                           control.method.list = list(name = "Stephenson", s = 6))
  
  m2_treat.s2_control.s2 = method_combine(dt$Z, dt$Y_obs, N = n, k_vec = k_vec, simul = simul,
                                           treat.method.list = list(name = "Stephenson", s = 2), 
                                           control.method.list = list(name = "Stephenson", s = 2))
  
  m2_treat.s6_control.s6 = method_combine(Z= dt$Z, Y = dt$Y_obs, N = n, k_vec = k_vec, simul = simul,
                                           treat.method.list = list(name = "Stephenson", s = 6), 
                                           control.method.list = list(name = "Stephenson", s = 6))
  
  m2_treat.s6_control.s2 = method_combine(Z= dt$Z, Y = dt$Y_obs, N = n, k_vec = k_vec, simul = simul,
                          treat.method.list = list(name = "Stephenson", s = 6), 
                          control.method.list = list(name = "Stephenson", s = 2))
  
  
  # Method M3 depends on simul
  m3_treat.s2_control.s6 = method_bergerboos(dt$Z, dt$Y_obs, N = n, k_vec = k_vec, simul = simul,
                                treat.method.list = list(name = "Stephenson", s = 2), 
                                control.method.list = list(name = "Stephenson", s = 6))
  m3_treat.s2_control.s2 = method_bergerboos(dt$Z, dt$Y_obs, N = n, k_vec = k_vec, simul = simul,
                            treat.method.list = list(name = "Stephenson", s = 2), 
                            control.method.list = list(name = "Stephenson", s = 2))
  m3_treat.s6_control.s6 = method_bergerboos(Z= dt$Z, Y = dt$Y_obs, N = n, k_vec = k_vec, simul = simul,
                            treat.method.list = list(name = "Stephenson", s = 6), 
                            control.method.list = list(name = "Stephenson", s = 6))
  m3_treat.s6_control.s2 = method_bergerboos(Z= dt$Z, Y = dt$Y_obs, N = n, k_vec = k_vec, simul = simul,
                             treat.method.list = list(name = "Stephenson", s = 6), 
                             control.method.list = list(name = "Stephenson", s = 2))
  
  
  res = rbind(m2_s2, m2_s6, 
    m3_treat.s2_control.s6, m3_treat.s2_control.s2, m3_treat.s6_control.s6, m3_treat.s6_control.s2, 
    m4_treat.s2_control.s6, m4_treat.s2_control.s2, m4_treat.s6_control.s6, m4_treat.s6_control.s2)
  
  res$method = c(rep('M1-S2', J), rep('M1-S6', J),
                 rep('M2-S2-S6', J), rep('M2-S2-S2', J), rep('M2-S6-S6', J),rep('M2-S6-S2', J),
                 rep('M3-S2-S6', J), rep('M3-S2-S2', J), rep('M3-S6-S6', J),rep('M3-S6-S2', J))

  return(res)
}

###############################################################################
dt_086 = read.csv("~/Chen_Zhe/HVTN086_data.csv")

run_simu <- function(nt, nc, antigen_of_interest  = 'Con 6 gp120/B', MC = 1000){
  n = nt+nc
  k_vec = c(ceiling(0.5*n), ceiling(0.75*n), 
            ceiling(0.8*n), ceiling(0.85*n),
            ceiling(0.90*n), ceiling(0.95*n))
  
  all_res_simul = NULL
  all_res_point = NULL
  
  for (i in 1:MC){
    dt = simulate_from_realdata(realdata = dt_086, antigen_of_interest = antigen_of_interest, 
                                arm.treated = "T1", arm.control ="T2", nt = nt, nc = nc)
    
    
    # Simultaneous Inference
    res_temp_simul = compare(dt, k_vec = k_vec, simul = TRUE)
    res_temp_simul$simultaneous = rep(TRUE, dim(res_temp_simul)[1])
    all_res_simul = rbind(all_res_simul, res_temp_simul)
    
    # Pointwise Inference
    res_temp_point = compare(dt, k_vec = k_vec, simul = FALSE)
    res_temp_point$simultaneous = rep(FALSE, dim(res_temp_point)[1])
    all_res_point = rbind(all_res_point, res_temp_point)
  }
  
  res = rbind(all_res_simul, all_res_point)
  res$nt = nt
  res$nc = nc
  res$antigen = antigen_of_interest
  res$tau = rep(sort(dt$tau)[k_vec], 2*length(unique(res$method)))
  return(res)
}

# RUN SIMULATIONS
res = NULL
for (nt in c(30, 50, 100)) {
    for (antigen in c('Con 6 gp120/B','gp41')){
        cat(nt, antigen, '\n')
        r = run_simu(nt=nt, nc=nt, antigen_of_interest  = antigen, MC = 1)
        res = rbind(res, r)
  }
}

# Cluster TASK_ID
touse<-as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
dir_name = paste0('/Results')
FileSave<-paste0(dir_name, '/sim_' , touse, ".csv")
write.table(res, file=FileSave, row.names = FALSE)

