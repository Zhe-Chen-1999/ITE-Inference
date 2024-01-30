setwd("~/Chen_Zhe/Real Data")
source("~/Chen_Zhe/functions.R")
data = read.csv('./dat_xprot_highResp.csv')

library(RIQITE)
library(dplyr)
library(ggplot2)
library(xtable)
# HVTN086
data_086 = data %>%
  filter(prot == 'v086') %>%
  filter(!antigen == 'ANY ENV')

data_086_treated_BAMA = data_086 %>%
  mutate(outcome = ifelse(outcome < 100, 100, outcome)) %>%
  dplyr::select(prot, rx_code, antigen, outcome, age, sex, race, bmi, 
                bpsys, bpdia, alt,
                hematocrit, hemoglobin, lymphocytes, 
                mcv, neutrophils, platelets, wbc_total, creatinine)


data_205 = data %>%
  filter(prot == 'v205') %>%
  filter(!antigen == 'ANY ENV') %>%
  mutate(outcome = ifelse(outcome < 100, 100, outcome)) %>%
  dplyr::select(prot, rx_code, antigen, outcome, age, sex, race, bmi,
                bpsys, bpdia, alt,
                hematocrit, hemoglobin, lymphocytes, 
                mcv, neutrophils, platelets, wbc_total, creatinine)


dt = rbind(data_086_treated_BAMA, data_205)


ggplot(data = dt %>%
         filter(antigen == 'gp41') %>%
         filter((prot == 'v086' & rx_code == 'T1') | (prot == 'v205' & rx_code == 'T4')),
       aes(y = outcome, x = prot)) + geom_boxplot() +
  geom_jitter(width = 0.2, shape = 2, color = 'steel blue', size = 1.5, stroke = 1.1) +
  facet_wrap(~antigen) + scale_y_continuous(limits = c(100, 35000), trans = 'log10',
                                            breaks = c(100, 1000, 2000, 5000, 10000, 20000, 35000)) +
  theme_bw(base_size = 25) + xlab('Treatment group') + ylab('Net MFI response')

dt = dt %>%
  mutate(outcome = log10(outcome))

########################## v205 T4 matched to 086 T1 ###############################################
library(match2C)
antigen_of_interest = 'gp41'

#table(dt$race) # all Black
dt$sex = ifelse(dt$sex == "Male", 1, 0)
data_analysis = dt %>%
  filter(antigen == antigen_of_interest) %>%
  filter((prot == 'v086' & rx_code == 'T1') | (prot == 'v205' & rx_code == 'T4')) %>%
  mutate(trt = ifelse(prot == 'v086' & rx_code == 'T1', 1, 0))

nt = sum(data_analysis$trt)
nc = sum(1-data_analysis$trt)
permute_id = sample(nt+nc, replace = FALSE)

table(data_analysis$race[which(data_analysis$prot == 'v205')])

X = data_analysis %>%
  dplyr::select(age, sex, 
                bmi, bpsys, bpdia, 
               hematocrit, hemoglobin, lymphocytes, 
               mcv, neutrophils, platelets) # covariates to be matched
Z = data_analysis$trt 

# Fit a propensity score model
propensity = glm(Z~age+
                   sex+bmi+bpsys + bpdia +
                   hematocrit + hemoglobin + lymphocytes + 
                 mcv + neutrophils + platelets,
                 family=binomial, data = X)$fitted.values
data_analysis$propensity = propensity

# Perform a matching with minimal input
matching_output = match_2C(Z = Z, X = X, 
                           propensity = propensity, controls = 1,
                           dataset = data_analysis)

# Matched dataset
cov_list = c('age', 
             'sex', 'bmi', 'bpsys', 'bpdia',
             'hematocrit', 'hemoglobin', 'lymphocytes', 
             'mcv', 'neutrophils', 'platelets', 'propensity')
matched_data = matching_output$matched_data_in_order %>%
  filter(!is.na(matched_set)) %>%
  dplyr::select(all_of(cov_list), matched_set, trt)

# Matched treated
dt_treated_after = matched_data %>% filter(trt == 1)

# Matched control
dt_control_after = matched_data %>% filter(trt == 0)


###### Balance table ######
# Before matching
d = data_analysis
dt_treated_before <- d[d$trt == 1, cov_list]
dt_control_before <- d[d$trt == 0, cov_list]

mean_treated_before = apply(dt_treated_before, 2, mean)
mean_control_before = apply(dt_control_before, 2, mean)
mean_diff_before = mean_treated_before - mean_control_before

sd_treated_before = apply(dt_treated_before, 2, stats::sd)
sd_control_before = apply(dt_control_before, 2, stats::sd)
pooled_sd = sqrt(sd_treated_before^2 + sd_control_before^2)
std_before = mean_diff_before/pooled_sd

# After matching
mean_treated_after = apply(dt_treated_after[, cov_list], 2, mean)
mean_control_after = apply(dt_control_after[, cov_list], 2, mean)
mean_diff_after = mean_treated_after - mean_control_after

std_after = mean_diff_after/pooled_sd
balance_table = data.frame(mean_treated_before, mean_control_before, 
                           std_before, mean_control_after, std_after)
print(xtable(balance_table, digits=2))
View(balance_table)

###### Plot propensity score ######
propens_treated = d$propensity[d$trt == 1]
propens_all_control = d$propensity[d$trt == 0]
propens_matched_control = dt_control_after$propensity
group = c(rep("Treated", length(propens_treated)), 
          rep("All controls", length(propens_all_control)), 
          rep("Matched controls", length(propens_matched_control)))
propensity = c(propens_treated, propens_all_control, 
               propens_matched_control)
propens_df = data.frame(group, propensity)
pt = ggplot2::ggplot(data = propens_df, 
                     ggplot2::aes(x = propensity, color = group)) + 
  ggplot2::geom_density(linewidth = 1.5) + 
  ggplot2::theme_bw(base_size = 20) + 
  ggplot2::theme(legend.position = "top")

print(pt)

########################## Stratified CRE: 086.T1 (treated) vs 205.T4 (control) ##########################
# devtools::install_github("Yongchang-Su/QIoT")

library(QIoT)
matched_dt = matching_output$matched_data_in_order %>%
  filter(!is.na(matched_set))
head(matched_dt, 6)

matched_dt$lab = ifelse(matched_dt$rx_code == 'T1', 'HVTN 086 T1', 'HVTN 205 T4')
boxplot_205T4_086T1 = ggplot(data = matched_dt,
       aes(y = 10^outcome, x = lab)) +
  geom_jitter(width = 0.2, shape = 2, color = 'steel blue', size = 1.5, stroke = 1.1) +
  facet_wrap(~antigen) + scale_y_continuous(limits = c(4500, 35000), trans = 'log10',
                                            breaks = c(4500, 6000, 10000, 20000, 35000)) +
  theme_bw(base_size = 30) + xlab('Treatment group') + ylab('Net MFI response')

ggsave('boxplot_205T4_086T1.pdf', plot = boxplot_205T4_086T1, 
       path =  "/home/bzhang3/Chen_Zhe/Real Data",
       scale = 1, width = 12, height = 8, 
       units = "in")

S = 33
n = 2
N = S*n
permute_id = permute(1:N)
block = matched_dt$matched_set[permute_id]
Z = matched_dt$trt[permute_id] # T1 086: Z=1; T4 205: Z=0
Y = matched_dt$outcome[permute_id]
boxplot(Y[Z==1],Y[Z==0])
wilcox.test(Y[Z==1],Y[Z==0],conf.int = TRUE, paired = TRUE, alternative = "greater")$conf.int[1] #95 percent confidence interval: [0.07, Inf) #sample estimates: (pseudo)median 0.09207666 

method.list.all = list()
method.list.all[[1]] = list(name = "Wilcoxon") # s=2 is same as method.list.all[[1]] = list(name = "Wilcoxon")

### Calculate lower limit of prediction intervals for treated units
CIs = ci_quantile_scre(Z, Y, block, method.list.all = method.list.all, 
                       alternative = "upper", confidence = 0.975, switch = FALSE)
CI_lower.treat = CIs$LB[34:66] 

### Calculate lower limit of prediction intervals for control units
CI.control = ci_quantile_scre(1-matched_dt$trt[permute_id], -matched_dt$outcome[permute_id], block, method.list.all = method.list.all, 
                       alternative = "upper", confidence = 0.975)
CI_lower.control = CI.control$LB[34:66] # 0.05662749 - 0.05662807

### Combine inference for two group
ci.all = sort( c(CI_lower.treat, CI_lower.control) )

pp = data.frame(pp_0 = sum(ci.all > 0), 
             pp_0.5 = sum(ci.all > 0.5), 
             pp_1 = sum(ci.all > 1), 
             pp_1.5 = sum(ci.all > 1.5))

ci_lower = data.frame(x = ci.all[ci.all>-Inf], y = (1:length(ci.all))[ci.all>-Inf])
xend = ceiling(max(ci_lower$x))

ci_plot_s2 = ggplot(data = ci_lower) + geom_segment(data = ci_lower, aes(x = x, y = y, xend = 0.1, yend = y))+
  geom_point(aes(x = x, y = y), size = 2.5) +
  scale_y_continuous(limits = c(43, N),  breaks = seq(44,N,2)) +
  theme_bw(base_size = 30) +
  xlab(expression("Magnitude of treatment effect"~tau[(k)])) + ylab('k')

ggsave('086T1_205T4_scre_Gamma1.pdf', plot = ci_plot_s2, 
       #path = '',
       scale = 1, width = 12, height = 8, 
       units = "in")


########################## Sensitivity Analysis: 086.T1 (treated) vs 205.T4 (control) ##########################
## CIs for all units under different Gamma
for (gamma in c(1.2, 1.5, 2.5, 3.3)) {

  CI_lower.treat = ci_quantile_sen(Z, Y, block, gam=gamma, quantiles= (nt+1):N,
                                 alternative = "upper", method.list.all=method.list.all, confidence = 0.975)$LB

  CI_lower.control = ci_quantile_sen(1-Z, -Y, block, gam=gamma, quantiles= (nt+1):N,
                                   alternative = "upper", method.list.all=method.list.all, confidence = 0.975)$LB

  CI_lower = sort( c(CI_lower.treat, CI_lower.control) )
  all_ci = data.frame(x = CI_lower, y = (1:length(CI_lower)), gam = rep(gamma,length(CI_lower)))
  saveRDS(all_ci, file = paste0("ci_lower_", gamma,".rds"))
}

## CI plots under different Gamma
for (gamma in c(1.2, 1.5, 2.5, 3.3)) {
  ci_lower = readRDS(file = paste0("ci_lower_", gamma,".rds"))
  ci_lower$x_new =  ifelse(is.infinite(ci_lower$x), -0.05, ci_lower$x)
  p = ggplot(data = ci_lower) + geom_segment(data = ci_lower, aes(x = x_new, y = y, xend = 0.1, yend = y))+
    geom_point(data = ci_lower[!is.infinite(ci_lower$x),],aes(x = x, y = y), size = 2.5) +
    scale_x_continuous(limits = c(-0.05, 0.1))+
    scale_y_continuous(limits = c(47, N),  breaks = seq(47+1,N,2)) +
    theme_bw(base_size = 30) +
    xlab(expression("Magnitude of treatment effect"~tau[(k)])) + ylab('k')
  
  ggsave(paste0("086T1_205T4_Gamma",gamma,".pdf"), plot = p, 
         path = "~/Chen_Zhe/Real Data/SA",
         scale = 1, width = 12, height = 8, 
         units = "in")
}

## proportion of units with positive effects under different Gamma
for (gamma in c(1.2, 1.5, 2.5, 3.3)) {
  ci_lower = readRDS(file = paste0("ci_lower_", gamma,".rds"))
  print(paste0(gamma," ", sum(ci_lower$x>0)/N))
}

ci_lower_gamma1.5 = readRDS(file = paste0("ci_lower_", 1.5,".rds"))
ci_lower_gamma1.5$x>0

ci_lower_gamma3.3 = readRDS(file = paste0("ci_lower_", 3.3,".rds"))
ci_lower_gamma3.3$x>0
