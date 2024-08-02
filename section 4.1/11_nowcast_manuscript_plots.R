# > sessionInfo()
# R version 4.4.1 (2024-06-14) -- "Race for Your Life"
# Platform: aarch64-apple-darwin20
# Running under: macOS Monterey 12.4

library(dplyr)     # dplyr_1.1.4
library(ggplot2)   # ggplot2_3.5.1
library(lubridate) #lubridate_1.9.3
library(cowplot)   #cowplot_1.1.3

############################
## Now cast first wave plot
############################

load(paste0('./nowcast_firstwave/results_allmodels.RData'))
load(paste0('./nowcast_firstwave/CANresults_allmodels.RData'))

load(paste0('./nowcast_firstwave_xiaotian/results_allmodels_xiaotian.RData'))
load(paste0('./nowcast_firstwave_xiaotian/CANresults_allmodels_xiaotian.RData'))


gg1 <- df_full_can_firstwave %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, ave_exp_v_u_fixed, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_errorbar(aes(ymax = ave_exp_v_u_fixed_upr, ymin = ave_exp_v_u_fixed_lwr), show.legend = FALSE, size= 0.5)+
  geom_ribbon(data =df_full_can_firstwave %>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
                filter(sample_date >= ymd("2021-11-20") & sample_date <= ymd("2022-01-19")),
              aes(x = sample_date, ymax = ave_exp_v_u_fixed_upr, ymin =ave_exp_v_u_fixed_lwr), alpha = 0.1, size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)")), breaks = scales::pretty_breaks(n=5),limits = c(0,800))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 18))


gg2 <- df_full_can_firstwave %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, ave_exp_v_u_fixed_deriv, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_errorbar(aes(ymax = ave_exp_v_u_fixed_deriv_upr, ymin = ave_exp_v_u_fixed_deriv_lwr), show.legend = FALSE, size= 0.5)+
  geom_ribbon(data =df_full_can_firstwave %>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
                filter(sample_date >= ymd("2021-11-20") & sample_date <= ymd("2022-01-19")),
              aes(x = sample_date, ymax = ave_exp_v_u_fixed_deriv_upr, ymin =ave_exp_v_u_fixed_deriv_lwr), alpha = 0.1, size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)'")), breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 18))


gg3 <- df_full_can_firstwave %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, post_prob_ave_exp_v_u_fixed_deriv, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_line(data =df_full_can_firstwave %>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
              filter(sample_date >= ymd("2021-11-20") & sample_date <= ymd("2022-01-19")),
            aes(x = sample_date, y=post_prob_ave_exp_v_u_fixed_deriv), size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = "Probability", breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 14))



gg4 <- df_full_can_firstwave_xiaotian %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, p1_exp_increase, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_line(data =df_full_can_firstwave_xiaotian %>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
              filter(sample_date >= ymd("2021-11-20") & sample_date <= ymd("2022-01-19")),
            aes(x = sample_date, y=p1_exp_increase), size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = "Probability", breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 14))


gg5 <- df_full_can_firstwave_xiaotian %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, ave_exp_f1_IS_fixed, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_errorbar(aes(ymax = ave_exp_f1_IS_fixed_upr, ymin = ave_exp_f1_IS_fixed_lwr), show.legend = FALSE, size= 0.5)+
  geom_ribbon(data =df_full_can_firstwave_xiaotian %>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
                filter(sample_date >= ymd("2021-11-20") & sample_date <= ymd("2022-01-19")),
              aes(x = sample_date, ymax = ave_exp_f1_IS_fixed_upr, ymin =ave_exp_f1_IS_fixed_lwr), alpha = 0.1, size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)")), breaks = scales::pretty_breaks(n=5),limits = c(0,800))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 18))




fest1 = cowplot::plot_grid(add_sub(gg1,"a) Proposed model, station average signal",size = 12),
                           add_sub(gg3,"c) Proposed model, probability of increase",size = 12),
                           add_sub(gg2,"e) Proposed model, derivative of station average signal",size = 12),
                           add_sub(gg5,"b) Comparison model, station average signal",size = 12),
                           add_sub(gg4,"d) Comparison model, probability of increase",size = 12),
                           ncol=2, align="v", byrow = FALSE) 




############################
## Now cast no wave plot
############################ 


load(paste0('./nowcast_nowave/results_allmodels.RData'))
load(paste0('./nowcast_nowave/CANresults_allmodels.RData'))

load(paste0('./nowcast_nowave_morebasis_xiaotian/results_allmodels_xiaotian.RData'))
load(paste0('./nowcast_nowave_morebasis_xiaotian/CANresults_allmodels_xiaotian.RData'))

min_date = "2022-12-01"
max_date = "2023-03-30"

gg1 <- df_full_can_nowave %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, ave_exp_v_u_fixed, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_errorbar(aes(ymax = ave_exp_v_u_fixed_upr, ymin = ave_exp_v_u_fixed_lwr), show.legend = FALSE, size= 0.5)+
  geom_ribbon(data =df_full_can_nowave %>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
                filter(sample_date >= min_date & sample_date <= max_date),
              aes(x = sample_date, ymax = ave_exp_v_u_fixed_upr, ymin =ave_exp_v_u_fixed_lwr), alpha = 0.1, size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)")), breaks = scales::pretty_breaks(n=5),limits = c(0,200))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 18))



gg2 <- df_full_can_nowave %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, ave_exp_v_u_fixed_deriv, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_errorbar(aes(ymax = ave_exp_v_u_fixed_deriv_upr, ymin = ave_exp_v_u_fixed_deriv_lwr), show.legend = FALSE, size= 0.5)+
  geom_ribbon(data =df_full_can_nowave %>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
                filter(sample_date >= min_date & sample_date <= max_date),
              aes(x = sample_date, ymax = ave_exp_v_u_fixed_deriv_upr, ymin =ave_exp_v_u_fixed_deriv_lwr), alpha = 0.1, size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)'")), breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 18))


gg3 <- df_full_can_nowave %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, post_prob_ave_exp_v_u_fixed_deriv, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_line(data =df_full_can_nowave%>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
              filter(sample_date >= min_date & sample_date <= max_date),
            aes(x = sample_date, y=post_prob_ave_exp_v_u_fixed_deriv), size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = "Probability", breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 14))


gg4 <- df_full_can_nowave_xiaotian %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, p1_exp_increase, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_line(data =df_full_can_nowave_xiaotian %>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
              filter(sample_date >= min_date & sample_date <= max_date),
            aes(x = sample_date, y=p1_exp_increase), size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = "Probability", breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 14))


gg5 <- df_full_can_nowave_xiaotian %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, ave_exp_f1_IS_fixed, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_errorbar(aes(ymax = ave_exp_f1_IS_fixed_upr, ymin = ave_exp_f1_IS_fixed_lwr), show.legend = FALSE, size= 0.5)+
  geom_ribbon(data =df_full_can_nowave_xiaotian %>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
                filter(sample_date >= min_date & sample_date <= max_date),
              aes(x = sample_date, ymax = ave_exp_f1_IS_fixed_upr, ymin =ave_exp_f1_IS_fixed_lwr), alpha = 0.1, size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)")), breaks = scales::pretty_breaks(n=5),limits = c(0,200))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18))


fest1 = cowplot::plot_grid(add_sub(gg1,"a) Proposed model, station average signal",size = 12),
                           add_sub(gg3,"c) Proposed model, probability of increase",size = 12),
                           add_sub(gg2,"e) Proposed model, derivative of station average signal",size = 12),
                           add_sub(gg5,"b) Comparison model, station average signal",size = 12),
                           add_sub(gg4,"d) Comparison model, probability of increase",size = 12),
                           ncol=2, align="v", byrow = FALSE) 



