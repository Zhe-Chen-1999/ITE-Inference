setwd("~/Chen_Zhe/Results/reviewpaper_simulation_from_realdata_1211")
all_res = read.csv('merged_res_for_1000_simulations.csv')
all_res$method = ifelse(all_res$method == 'M2-S10', 'M2-S6', all_res$method)

################### Median Plot##########################
median_tab <- all_res %>%
  group_by(antigen, simultaneous, nt, method, k) %>%
  summarise(median = median(lower)) 

# median_table <- median_tab %>% pivot_wider(names_from = method, values_from = median)
# View(median_table)
# summary(median_tab$median[(median_tab$nt==30) & (!is.infinite(median_tab$median))])

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
    mutate(method = ifelse(method =="M2-S2", "M1-S2", 
                           ifelse(method =="M2-S6", "M1-S6",
                                  ifelse(method =="M3-t.S2-c.S6", "M2-S2-S6", 
                                         ifelse(method =="M3-t.S2-c.S2", "M2-S2-S2", 
                                                ifelse(method =="M3-t.S6-c.S6", "M2-S6-S6", 
                                                       ifelse(method =="M3-t.S6-c.S2", "M2-S6-S2", 
                                                              ifelse(method =="M4-t.S2-c.S6", "M3-S2-S6", 
                                                                     ifelse(method =="M4-t.S2-c.S2", "M3-S2-S2", 
                                                                            ifelse(method =="M4-t.S6-c.S6", "M3-S6-S6", 
                                                                                   ifelse(method =="M4-t.S6-c.S2", "M3-S6-S2", method)))))))))))%>%
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
         path = '/home/bzhang3/Chen_Zhe/Results/reviewpaper_simulation_from_realdata_1211/update_0129',
         scale = 1, width = 12, height = 10, 
         units = "in")
}


####################### measurement II: SS #####################################
all_res$lower[is.infinite(all_res$lower)] = -10

SS_avg <- all_res %>%
  #filter(simultaneous == TRUE) %>%
  mutate(nc = nt)%>%
  group_by(antigen, nt, nc, simultaneous, method) %>%
  summarise(SS_avg = round(mean((lower - tau)^2), digits=3)) 

library(tidyr)
SS_table <- SS_avg %>% pivot_wider(names_from = method, values_from = SS_avg)
View(SS_table)

library(xtable)
tb <- SS_table%>%
  filter(p == 0.8)
tb <- tb[,-c(1,2)]
tb <- tb[,c(1,2,3,7,6,4,5,11,10,8,9)]
xtable(SS_table)

print(xtable(SS_table), include.rownames=FALSE)
