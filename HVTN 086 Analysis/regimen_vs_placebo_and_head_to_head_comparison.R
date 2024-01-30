source("~/helper_functions.R")

library(RIQITE)
library(dplyr)
library(ggplot2)

#################### Processing data ################################
# HVTN086
data = read.csv('./dat_xprot_highResp.csv')
data_086 = data %>%
  filter(prot == 'v086') %>%
  filter(!antigen == 'ANY ENV')

# Set BAMA limit of quantification
data_086_treated_BAMA = data_086 %>%
  mutate(outcome = ifelse(outcome < 100, 100, outcome)) %>%
  dplyr::select(rx_code, antigen, outcome)

# HVTN086 enrolled eight controls to accompany each treatment arm
# The outcomes are set to 100
data_086_control_BAMA = data.frame(rx_code = rep(c('C1', 'C2', 'C3', 'C4'), each = 24),
                                   antigen = rep(c('Con 6 gp120/B', 'Con S gp140 CFI', 'gp41'), 8),
                                   outcome = 100)

ggplot(data = data_086_treated_BAMA, aes(y = outcome, x = rx_code)) + geom_boxplot() +
  geom_jitter(width = 0.2, shape = 2, color = 'steel blue', size = 1.5, stroke = 1.1) +
  facet_wrap(~antigen) + scale_y_continuous(limits = c(100, 35000), trans = 'log10',
                                            breaks = c(100, 1000, 2000, 5000, 10000, 20000, 35000)) +
  theme_bw(base_size = 25) + xlab('Treatment group') + ylab('Net MFI response')


data_086_BAMA = bind_rows(data_086_treated_BAMA, data_086_control_BAMA)

dt = data_086_BAMA %>%
  mutate(outcome = log10(outcome))

################### Additional helper functions ##################################

CL_construct_all_CI <- function(dt, antigen_of_interest, arm_1, arm_2, method, simul = FALSE, 
                                caughey.method.list = list(name = "Stephenson", s = 2),
                                treat.method.list= list(name = "Stephenson", s = 2), control.method.list = list(name = "Stephenson", s = 10)){
  data_analysis = dt %>%
    filter(antigen == antigen_of_interest) %>%
    filter(rx_code == arm_1 | rx_code == arm_2) %>%
    mutate(trt = ifelse(rx_code == arm_1, 1, 0))
  nt = sum(data_analysis$trt)
  nc = sum(1-data_analysis$trt)
  permute_id = sample(nt+nc, replace = FALSE)
  
  all_ci = NULL
  if(method == 'S&M'){
    all_ci = method_SM(Z = data_analysis$trt[permute_id], 
                       Y = data_analysis$outcome[permute_id], 
                       LOD = rep(2,nt+nc), k_vec = 1:(nt+nc), simul = simul)
  }else if(method == 'M1'){
    all_ci = method_caughey(data_analysis$trt[permute_id], data_analysis$outcome[permute_id], k_vec = 1:(nt+nc), method.list = caughey.method.list)
  }else if(method == 'M2'){
    all_ci =  method_combine(Z= data_analysis$trt[permute_id],
                             Y = data_analysis$outcome[permute_id],
                             N = nt+nc, 
                             treat.method.list = treat.method.list, control.method.list = control.method.list, 
                             k_vec = 1:(nt+nc), 
                             simul = simul)
    
  }else if(method == 'M3'){
    all_ci = method_bergerboos(Z = data_analysis$trt[permute_id], 
                      Y = data_analysis$outcome[permute_id],
                      N = nt+nc, 
                      k_vec = 1:(nt+nc), 
                      treat.method.list = treat.method.list, control.method.list = control.method.list, 
                      simul = simul)
  }
  return(all_ci)
}

visualize_CI <- function(all_ci, mean_diff_lb, q_vec, xlim_L = -2.5, xlim_R = 3, xend = 3){
  
  plot_ci = all_ci %>%
    group_by(gp) %>%
    filter(k %in% ceiling(n() * q_vec)) %>%
    mutate(k = q_vec)
  
  plot_df = data.frame(x = plot_ci$lower, y = plot_ci$k,
                       xend = xend, yend = plot_ci$k,
                       gp = plot_ci$gp)
  
  ci_plot = ggplot(data = plot_df) + geom_segment(data = plot_df, aes(x = x, y = y, xend = xend, yend = yend))+
    geom_point(aes(x = x, y = y), size = 2.5) + xlim(c(xlim_L, xlim_R)) +
    geom_vline(data = mean_diff_lb, aes(xintercept =  diff), size = 1.5, linetype = 'dashed', color = 'red', alpha = 0.5) +
    scale_y_continuous(labels = scales::percent) + theme_bw(base_size = 25) +
    facet_wrap(~gp, ncol = 2)+
    xlab('Magnitude of treatment effect') + ylab('Quantile')
  
  return(ci_plot)
}

###########################################################################################
############  Characterizing immunogenicity profiles against placebo ######################
###########################################################################################

antigen_of_interest = 'Con 6 gp120/B'

############ pointwise inference ############ 
ci_1 = CL_construct_all_CI(dt, antigen_of_interest, arm_1 = 'T1', arm_2 = 'C1', method = 'S&M', simul = FALSE) %>%
  mutate(gp = 'T1 (n = 33) versus C1 (n = 8)')
ci_2 = CL_construct_all_CI(dt, antigen_of_interest, arm_1 = 'T2', arm_2 = 'C2', method = 'S&M', simul = FALSE) %>%
  mutate(gp = 'T2 (n = 32) versus C2 (n = 8)')
ci_3 = CL_construct_all_CI(dt, antigen_of_interest, arm_1 = 'T3', arm_2 = 'C3', method = 'S&M', simul = FALSE) %>%
  mutate(gp = 'T3 (n = 24) versus C3 (n = 8)')
ci_4 = CL_construct_all_CI(dt, antigen_of_interest, arm_1 = 'T4', arm_2 = 'C4', method = 'S&M', simul = FALSE) %>%
  mutate(gp = 'T4 (n = 23) versus C4 (n = 8)')
ci = bind_rows(ci_1, ci_2, ci_3, ci_4)

## pointwise lower confidence limit for N(c), with c = 0, 0.5, 1, 1.5, 2
pp <- ci %>% group_by(gp) %>%
  summarise(pp_0 = sum(lower > 0), pp_0.5 = sum(lower > 0.5), pp_1 = sum(lower > 1), pp_1.5 = sum(lower > 1.5), pp_2 = sum(lower > 2))

print(xtable(pp), include.rownames=FALSE)

## Confidence Interval For Difference of Means: one sample t-test
s = dt %>% filter(antigen == antigen_of_interest)
x1 = s$outcome[which(s$rx_code=='T1')]
x2 = s$outcome[which(s$rx_code=='T2')]
x3 = s$outcome[which(s$rx_code=='T3')]
x4 = s$outcome[which(s$rx_code=='T4')]

T1meanLB = t.test(x1, alternative = 'greater', conf.level = 0.95)$conf.int[1]
T2meanLB = t.test(x2, alternative = 'greater', conf.level = 0.95)$conf.int[1]
T3meanLB = t.test(x3, alternative = 'greater', conf.level = 0.95)$conf.int[1]
T4meanLB = t.test(x4, alternative = 'greater', conf.level = 0.95)$conf.int[1]

mean_diff_lb = tibble(gp = c('T1 (n = 33) versus C1 (n = 8)',
                             'T2 (n = 32) versus C2 (n = 8)',
                             'T3 (n = 24) versus C3 (n = 8)',
                             'T4 (n = 23) versus C4 (n = 8)'),
                      diff = c(T1meanLB-2, T2meanLB-2, T3meanLB-2, T4meanLB-2))

## CI plot
pointwise_ci_plot = visualize_CI(ci, mean_diff_lb, 
                  q_vec = c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95),
                  xlim_L = 0, xlim_R = 2.5, xend = 2.5)+
  ggtitle('IgG Responses to Con 6 gp120/B') +
  theme(plot.title=element_text(hjust=0.5))

pointwise_ci_plot

ggsave('Con6_regimen_vs_placebo_pointwise.pdf', plot = pointwise_ci_plot, 
       path = "~/Chen_Zhe/SCID_real_data",
       scale = 1, width = 12, height = 8, 
       units = "in")

### the largest ITE is at least 2.49
ci_1[ci_1$k == nrow(ci_1), "lower"]
### the top 25% ITE is at least 2.40
ci_1[ci_1$k == ceiling(0.75*nrow(ci_1)), "lower"]
### the median ITE is at least 2.30
ci_1[ci_1$k == ceiling(0.5*nrow(ci_1)), "lower"]

############ simultaneous inference ############ 
ci_1 = CL_construct_all_CI(dt, antigen_of_interest, arm_1 = 'T1', arm_2 = 'C1', method = 'S&M', simul = TRUE) %>%
  mutate(gp = 'T1 (n = 33) versus C1 (n = 8)')
ci_2 = CL_construct_all_CI(dt, antigen_of_interest, arm_1 = 'T2', arm_2 = 'C2', method = 'S&M', simul = TRUE) %>%
  mutate(gp = 'T2 (n = 32) versus C2 (n = 8)')
ci_3 = CL_construct_all_CI(dt, antigen_of_interest, arm_1 = 'T3', arm_2 = 'C3', method = 'S&M', simul = TRUE) %>%
  mutate(gp = 'T3 (n = 24) versus C3 (n = 8)')
ci_4 = CL_construct_all_CI(dt, antigen_of_interest, arm_1 = 'T4', arm_2 = 'C4', method = 'S&M', simul = TRUE) %>%
  mutate(gp = 'T4 (n = 23) versus C4 (n = 8)')
ci = bind_rows(ci_1, ci_2, ci_3, ci_4)
  
simul_ci_plot = visualize_CI(ci, mean_diff_lb, 
                  q_vec = c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95),
                  xlim_L = 0, xlim_R = 2.5, xend = 2.5)+
  ggtitle('IgG Responses to Con 6 gp120/B') +
  theme(plot.title=element_text(hjust=0.5))
simul_ci_plot

ggsave('Con6_regimen_vs_placebo_simul.pdf', plot = simul_ci_plot, 
       path = "~/Chen_Zhe/SCID_real_data",
       scale = 1, width = 12, height = 8, 
       units = "in")

###########################################################################################
###################### Head-to-head comparisons of two vaccine regimens ###################
###########################################################################################

###### simultaneous inference for gp41 using M2-S2-S6 ######

antigen_of_interest = 'gp41'

ci_gp41_T1_vs_T2 = CL_construct_all_CI(dt, antigen_of_interest, arm_1 = 'T1', arm_2 = 'T2', method = 'M2', simul = TRUE, 
                           treat.method.list = list(name = "Stephenson", s = 2), control.method.list = list(name = "Stephenson", s = 6)) %>%
  mutate(gp = 'T1 (n = 33) versus T2 (n = 32)')
ci_gp41_T1_vs_T3 = CL_construct_all_CI(dt, antigen_of_interest, arm_1 = 'T1', arm_2 = 'T3', method = 'M2', simul = TRUE, 
                           treat.method.list = list(name = "Stephenson", s = 2), control.method.list = list(name = "Stephenson", s = 6)) %>%
  mutate(gp = 'T1 (n = 33) versus T3 (n = 24)')
ci_gp41_T1_vs_T4 = CL_construct_all_CI(dt, antigen_of_interest, arm_1 = 'T1', arm_2 = 'T4', method = 'M2', simul = TRUE, 
                           treat.method.list = list(name = "Stephenson", s = 2), control.method.list = list(name = "Stephenson", s = 6)) %>%
  mutate(gp = 'T1 (n = 33) versus T4 (n = 23)')
simul_ci_compare_two_regimen_gp41 = bind_rows(ci_gp41_T1_vs_T2, ci_gp41_T1_vs_T3, ci_gp41_T1_vs_T4)

###### 95% lower confidence limits of the mean difference derived from the two-sample Wilcoxon rank sum test ###### 
s = dt %>% filter(antigen == antigen_of_interest)
x = s$outcome[which(s$rx_code=='T1')]
y2 = s$outcome[which(s$rx_code=='T2')]
y3 = s$outcome[which(s$rx_code=='T3')]
y4 = s$outcome[which(s$rx_code=='T4')]

T1Mean_substract_T2Mean_LB = wilcox.test(x, y2, conf.int = TRUE, alternative = "greater")$conf.int[1]
T1Mean_substract_T3Mean_LB = wilcox.test(x, y3, conf.int = TRUE, alternative = "greater")$conf.int[1]
T1Mean_substract_T4Mean_LB = wilcox.test(x, y4, conf.int = TRUE, alternative = "greater")$conf.int[1]

regimen_mean_diff_lb = tibble(gp = c('T1 (n = 33) versus T2 (n = 32)',
                                     'T1 (n = 33) versus T3 (n = 24)',
                                     'T1 (n = 33) versus T4 (n = 23)'),
                              diff = c(T1Mean_substract_T2Mean_LB, T1Mean_substract_T3Mean_LB, T1Mean_substract_T4Mean_LB))

######## CI plot for gp41 ##############
CI_plot_gp41 = visualize_CI(simul_ci_compare_two_regimen_gp41, regimen_mean_diff_lb, 
                  q_vec = c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95),
                  xlim_L = -1, xlim_R = 2, xend = 2)+
                  ggtitle('IgG Responses to gp41') +
                  theme(plot.title=element_text(hjust=0.5))
CI_plot_gp41

ggsave('~/Chen_Zhe/Real Data/head_to_head_comparison/M3_t.S2_c.S6_gp41_head_to_head_comparison.pdf', plot = CI_plot_gp41, 
       scale = 1, width = 12, height = 8, 
       units = "in")

## Proportion of participants benefited more from T1 compared to T2 by at least 0.57, 1 or 1.25 in the log10-scale
sum(ci_gp41_T1_vs_T2$lower>=T1Mean_substract_T2Mean_LB)/nrow(ci_gp41_T1_vs_T2)
sum(ci_gp41_T1_vs_T2$lower>=1)/nrow(ci_gp41_T1_vs_T2)
sum(ci_gp41_T1_vs_T2$lower>=1.25)/nrow(ci_gp41_T1_vs_T2)
###########################################################################################

###### simultaneous inference for Con 6 gp120/B using M2-S6-S6 ######
antigen_of_interest = 'Con 6 gp120/B'

ci_con6_T1_vs_T2 = CL_construct_all_CI(dt, antigen_of_interest, arm_1 = 'T1', arm_2 = 'T2', method = 'M2', simul = TRUE, 
                           treat.method.list = list(name = "Stephenson", s = 6), control.method.list = list(name = "Stephenson", s = 6)) %>%
  mutate(gp = 'T1 (n = 33) versus T2 (n = 32)')
ci_con6_T1_vs_T3 = CL_construct_all_CI(dt, antigen_of_interest, arm_1 = 'T1', arm_2 = 'T3', method = 'M2', simul = TRUE, 
                           treat.method.list = list(name = "Stephenson", s = 6), control.method.list = list(name = "Stephenson", s = 6)) %>%
  mutate(gp = 'T1 (n = 33) versus T3 (n = 24)')
ci_con6_T1_vs_T4 = CL_construct_all_CI(dt, antigen_of_interest, arm_1 = 'T1', arm_2 = 'T4', method = 'M2', simul = TRUE, 
                           treat.method.list = list(name = "Stephenson", s = 6), control.method.list = list(name = "Stephenson", s = 6)) %>%
  mutate(gp = 'T1 (n = 33) versus T4 (n = 23)')
simul_ci_compare_two_regimen_con6 = bind_rows(ci_con6_T1_vs_T2, ci_con6_T1_vs_T3, ci_con6_T1_vs_T4)

###### 95% lower confidence limits of the mean difference derived from the two-sample Wilcoxon rank sum test ###### 
s = dt %>% filter(antigen == antigen_of_interest)
x = s$outcome[which(s$rx_code=='T1')]
y2 = s$outcome[which(s$rx_code=='T2')]
y3 = s$outcome[which(s$rx_code=='T3')]
y4 = s$outcome[which(s$rx_code=='T4')]

T1Mean_substract_T2Mean_LB = wilcox.test(x, y2, conf.int = TRUE, alternative = "greater")$conf.int[1]
T1Mean_substract_T3Mean_LB = wilcox.test(x, y3, conf.int = TRUE, alternative = "greater")$conf.int[1]
T1Mean_substract_T4Mean_LB = wilcox.test(x, y4, conf.int = TRUE, alternative = "greater")$conf.int[1]

regimen_mean_diff_lb = tibble(gp = c('T1 (n = 33) versus T2 (n = 32)',
                                     'T1 (n = 33) versus T3 (n = 24)',
                                     'T1 (n = 33) versus T4 (n = 23)'),
                              diff = c(T1Mean_substract_T2Mean_LB, T1Mean_substract_T3Mean_LB, T1Mean_substract_T4Mean_LB))

######## CI plot for Con6 ##############
CI_plot_Con6 = visualize_CI(simul_ci_compare_two_regimen_con6, regimen_mean_diff_lb, 
                            q_vec = c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95),
                            xlim_L = -1, xlim_R = 2, xend = 2)

CI_plot_Con6 = CI_plot_Con6 + ggtitle('IgG Responses to Con 6 gp120/B') +
  theme(plot.title=element_text(hjust=0.5))

CI_plot_Con6

ggsave('~/Chen_Zhe/Real Data/head_to_head_comparison/M3_t.S6_c.S6_Con6_head_to_head_comparison.pdf', plot = CI_plot_Con6, 
       scale = 1, width = 12, height = 8, 
       units = "in")

## Proportion of participants benefited more from T1 compared to T2 by at least 1.61 in the log10-scale
sum(ci_con6_T1_vs_T2$lower>=T1Mean_substract_T2Mean_LB)/nrow(ci_con6_T1_vs_T2)

###########################################################################################
