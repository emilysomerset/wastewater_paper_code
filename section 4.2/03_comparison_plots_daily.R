# > sessionInfo()
# R version 4.4.1 (2024-06-14) -- "Race for Your Life"
# Platform: aarch64-apple-darwin20
# Running under: macOS Monterey 12.4

library(ggplot2)   #ggplot2_3.5.1
library(dplyr)     #dplyr_1.1.4 
library(lubridate) #lubridate_1.9.3
library(cowplot)   #cowplot_1.1.3

load(file = "results_comparisonmodel.RData")

gg1 = full_results %>% 
  ungroup() %>% 
  group_by(time_index) %>% 
  slice(1) %>% 
  ggplot(aes(Date,log_e_gc))+ 
  # geom_point(size = 0.5)+
  # facet_wrap(~uwwName, scales = "free")+ 
  geom_line(aes(Date, f1_int), col = "black",size = 0.5)+ 
  geom_ribbon(aes(Date, ymax = f1_int_upr, ymin = f1_int_lwr), alpha = 0.2, size = 0.5, fill = "black")+
  theme_bw()+
  scale_y_continuous(name = expression(v[t]), breaks = scales::pretty_breaks(n=5),
                     limits = c(6,12))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x),
                                   breaks = c(ymd("2021-09-01","2022-01-01"))))+ 
  theme(axis.title.y = element_text(size = 12))+
  theme(axis.text.x.top = element_text(vjust = -80),
        axis.ticks.x.top = element_blank() )
  
gg2 = full_results %>% 
  ungroup() %>% 
  group_by(time_index) %>% 
  slice(1) %>% 
  ggplot(aes(Date,log_e_gc))+ 
  # geom_point(size = 0.5)+
  # facet_wrap(~uwwName, scales = "free")+ 
  geom_line(aes(Date, f1_int_deriv, col = "black"),size = 0.5)+ 
  geom_ribbon(aes(Date, ymax = f1_int_deriv_upr, ymin = f1_int_deriv_lwr), alpha = 0.4, size = 0.5,col = "black", linetype = "dotted", fill = NA)+
  geom_line(aes(Date, f1_int_deriv_7/7, col = "red"),size = 0.5)+ 
  geom_ribbon(aes(Date, ymax = f1_int_deriv_7_upr/7, ymin = f1_int_deriv_7_lwr/7), alpha = 0.2, size = 0.5, fill = "red")+
  theme_bw()+
  scale_y_continuous(name = "Rate of Change", breaks = scales::pretty_breaks(n=5),
                     limits = c(-0.6,0.6))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x),
                                   breaks = c(ymd("2021-09-01","2022-01-01"))))+ 
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_colour_manual(name = "", values = c("black"="black","red"="red"), labels = c("Lag-1","Lag-7"))+
  theme(legend.position = c(0.99,0.99),
        legend.justification = c("right","top"),
        legend.margin = margin(-10,6,0,6))+
  theme(axis.title.y = element_text(size = 12))+
  theme(axis.text.x.top = element_text(vjust = -80),
        axis.ticks.x.top = element_blank())

load("results_ourmodel.RData")

df_full <- results$df_full

gg5 <- df_full %>% 
  group_by(sample_date) %>% 
  slice(1) %>%
  ggplot(aes(sample_date, v))+
  geom_line(col = "black",size = 0.5)+
  geom_ribbon(aes(ymin =  v_lwr, ymax =  v_upr), alpha = 0.2, size = 0.5, fill = "black")+
  theme_bw()+
  scale_y_continuous(name = "V(t)", breaks = scales::pretty_breaks(n=5),
                     limits = c(6,12))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x),
                                   breaks = c(ymd("2021-09-01","2022-01-01"))))+ 
  theme(axis.text.x.top = element_text(vjust = -80),
        axis.ticks.x.top = element_blank())

gg6 <- df_full %>% 
  group_by(sample_date) %>% 
  slice(1) %>%
  ggplot(aes(sample_date, v_deriv))+
  geom_line(col = "black",size = 0.5)+
  geom_ribbon(aes(ymin =   v_deriv_lwr, ymax =   v_deriv_upr), alpha = 0.2, size = 0.5, fill = "black")+
  theme_bw()+
  scale_y_continuous(name = "V(t)'", breaks = scales::pretty_breaks(n=5),
                     limits = c(-0.6,0.6))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x),
                                   breaks = c(ymd("2021-09-01","2022-01-01"))))+ 
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme(axis.text.x.top = element_text(vjust = -80),
        axis.ticks.x.top = element_blank())


fest1 = plot_grid(add_sub(gg5,"a) Proposed model, common signal",size = 10),
                  add_sub(gg1,"b) Comparison model, common signal",size = 10),
                  add_sub(gg6,"c) Proposed model, derivative of common signal",size = 10),
                  add_sub(gg2,"d) Comparison model, rates of change of common signal",size = 10))



