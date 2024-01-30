library(tidyr)

all_res = read.csv('merged_res_for_1000_simulations.csv')

################### Median Plot##########################
median_tab <- all_res %>%
  group_by(antigen, simultaneous, nt, method, k) %>%
  summarise(median = median(lower)) 

for (nt in c(30, 50 ,100)) {
  n = 2*nt
  k_vec = c(ceiling(0.95*n), ceiling(0.9*n), 
            ceiling(0.85*n), ceiling(0.8*n),
            ceiling(0.75*n), ceiling(0.5*n))
  
  median_plot = median_tab %>%
    mutate(median = ifelse(median == -Inf, -0.5, median)) %>% 
    mutate(simultaneous = ifelse(simultaneous == FALSE, "Pointwise Inference", "Simultaneous Inference")) %>% 
    filter(k %in% k_vec) %>%
    mutate(Quantile = ifelse(k == k_vec[1], '95 Percentile', 
                             ifelse(k == k_vec[2], '90 Percentile', 
                                    ifelse(k == k_vec[3], '85 Percentile', 
                                           ifelse(k == k_vec[4], '80 Percentile', 
                                                  ifelse(k == k_vec[5], '75 Percentile', ifelse(k == k_vec[6], '50 Percentile', k))))))) %>%
    ggplot(aes(x=method, y = median, group=Quantile, ymin=-0.5)) +
    geom_line(aes(color=Quantile))+ 
    geom_point(aes(color=Quantile))+
    facet_grid(antigen~simultaneous#,  scales="free_y"
               ) + 
    geom_point(data = . %>% group_by(antigen, simultaneous,Quantile)%>% filter(median == max(median)), shape=24, color="red",fill="red",size=3)+
    theme_bw(base_size = 22)+
    scale_color_viridis_d()+
    scale_y_continuous(breaks = c(-0.5, -0.35, 0, 0.5, 1, 1.5),labels=c(-Inf,-0.35, 0, 0.5, 1, 1.5))+
    xlab("Method")+
    ylab("Median of 95% lower confidence limits")+ 
    theme(legend.position = "top", panel.border = element_rect(colour = "black", fill=NA), axis.text.x = element_text(angle = 90))
  
  FileSave<-paste0("median_plot_nt", nt,'_nc', nt, ".pdf")
  ggsave(FileSave, plot = median_plot, 
         scale = 1, width = 12, height = 10, 
         units = "in")
}

####################### measurement II: Average SS #####################################

all_res$lower[is.infinite(all_res$lower)] = -10

SS_avg <- all_res %>%
  mutate(nc = nt)%>%
  group_by(antigen, nt, nc, simultaneous, method) %>%
  summarise(SS_avg = round(mean((lower - tau)^2), digits=3)) 

SS_table <- SS_avg %>% pivot_wider(names_from = method, values_from = SS_avg)
print(xtable(SS_table), include.rownames=FALSE)
