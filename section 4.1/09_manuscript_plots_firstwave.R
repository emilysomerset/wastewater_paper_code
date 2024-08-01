
load(paste0('./Xiaotian et al/phac/nowcast_firstwave2/results_allmodels.RData'))
load(paste0('./Xiaotian et al/phac/nowcast_firstwave2/CANresults_allmodels.RData'))

load(paste0('./Xiaotian et al/phac/nowcast_firstwave_xiaotian/results_allmodels_xiaotian.RData'))
load(paste0('./Xiaotian et al/phac/nowcast_firstwave_xiaotian/CANresults_allmodels_xiaotian.RData'))

## Rename some stuff for Xiaotian
# df_full_can_xiaotian <- df_full_can_xiaotian %>% 
#   rename("post_prob_ave_exp_f1_IS_fixed_deriv" = p1_exp_corr_increase)



### Get end point predictions. 

gg1 <- df_full_can_firstwave %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, ave_exp_f1_IS_fixed, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_errorbar(aes(ymax = ave_exp_f1_IS_fixed_upr, ymin = ave_exp_f1_IS_fixed_lwr), show.legend = FALSE, size= 0.5)+
  geom_ribbon(data =df_full_can_firstwave %>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
                filter(sample_date >= ymd("2021-11-20") & sample_date <= ymd("2022-01-19")),
              aes(x = sample_date, ymax = ave_exp_f1_IS_fixed_upr, ymin =ave_exp_f1_IS_fixed_lwr), alpha = 0.1, size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)")), breaks = scales::pretty_breaks(n=5),limits = c(0,800))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 18))

ggsave(filename = "./Xiaotian et al/phac/formanuscript_nowcast_firstwave/mubar_ourmodel.jpeg",
       plot = gg1, 
       device = "jpeg",
       width = 5, 
       height = 3,
       dpi = 300)



gg2 <- df_full_can_firstwave %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, ave_exp_f1_IS_fixed_deriv, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_errorbar(aes(ymax = ave_exp_f1_IS_fixed_deriv_upr, ymin = ave_exp_f1_IS_fixed_deriv_lwr), show.legend = FALSE, size= 0.5)+
  geom_ribbon(data =df_full_can_firstwave %>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
                filter(sample_date >= ymd("2021-11-20") & sample_date <= ymd("2022-01-19")),
              aes(x = sample_date, ymax = ave_exp_f1_IS_fixed_deriv_upr, ymin =ave_exp_f1_IS_fixed_deriv_lwr), alpha = 0.1, size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)'")), breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 18))

ggsave(filename = "./Xiaotian et al/phac/formanuscript_nowcast_firstwave/mubar_prime_ourmodel.jpeg",
       plot = gg2, 
       device = "jpeg",
       width = 5, 
       height = 3,
       dpi = 300)

gg3 <- df_full_can_firstwave %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, post_prob_ave_exp_f1_IS_fixed_deriv, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_line(data =df_full_can_firstwave %>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
              filter(sample_date >= ymd("2021-11-20") & sample_date <= ymd("2022-01-19")),
            aes(x = sample_date, y=post_prob_ave_exp_f1_IS_fixed_deriv), size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = "Probability", breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 14))

ggsave(filename = "./Xiaotian et al/phac/formanuscript_nowcast_firstwave/prob_ourmodel.jpeg",
       plot = gg3, 
       device = "jpeg",
       width = 5, 
       height = 3,
       dpi = 300)


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

ggsave(filename = "./Xiaotian et al/phac/formanuscript_nowcast_firstwave/prob_comparison.jpeg",
       plot = gg4, 
       device = "jpeg",
       width = 5, 
       height = 3,
       dpi = 300)


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

ggsave(filename = "./Xiaotian et al/phac/formanuscript_nowcast_firstwave/mubar_comparison.jpeg",
       plot = gg5, 
       device = "jpeg",
       width = 5, 
       height = 3,
       dpi = 300)

library(cowplot)
fest1 = cowplot::plot_grid(add_sub(gg1,"a) Proposed model, station average signal",size = 12),
                           add_sub(gg3,"c) Proposed model, probability of increase",size = 12),
                           add_sub(gg2,"e) Proposed model, derivative of station average signal",size = 12),
                           add_sub(gg5,"b) Comparison model, station average signal",size = 12),
                           add_sub(gg4,"d) Comparison model, probability of increase",size = 12),
                           ncol=2, align="v", byrow = FALSE) 


ggsave(filename = paste0("./Xiaotian et al/phac/formanuscript_nowcast_firstwave/comparebothmodels.jpeg"),
       plot = grid.arrange(fest1), 
       device = "jpeg",
       width = 8, 
       height = 8,
       dpi = 300)



## Add the other plot now. 


load(paste0('./Xiaotian et al/phac/nowcast_nowave/results_allmodels.RData'))
load(paste0('./Xiaotian et al/phac/nowcast_nowave/CANresults_allmodels.RData'))

load(paste0('./Xiaotian et al/phac/nowcast_nowave_morebasis_xiaotian/results_allmodels_xiaotian.RData'))
load(paste0('./Xiaotian et al/phac/nowcast_nowave_morebasis_xiaotian/CANresults_allmodels_xiaotian.RData'))

min_date = "2022-12-01"
max_date = "2023-03-30"

gg1 <- df_full_can %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, ave_exp_f1_IS_fixed, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_errorbar(aes(ymax = ave_exp_f1_IS_fixed_upr, ymin = ave_exp_f1_IS_fixed_lwr), show.legend = FALSE, size= 0.5)+
  geom_ribbon(data =df_full_can %>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
                filter(sample_date >= min_date & sample_date <= max_date),
              aes(x = sample_date, ymax = ave_exp_f1_IS_fixed_upr, ymin =ave_exp_f1_IS_fixed_lwr), alpha = 0.1, size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)")), breaks = scales::pretty_breaks(n=5),limits = c(0,200))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 18))

ggsave(filename = "./Xiaotian et al/phac/formanuscript_nowcast_nowave/mubar_ourmodel.jpeg",
       plot = gg1, 
       device = "jpeg",
       width = 5, 
       height = 3,
       dpi = 300)



gg2 <- df_full_can %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, ave_exp_f1_IS_fixed_deriv, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_errorbar(aes(ymax = ave_exp_f1_IS_fixed_deriv_upr, ymin = ave_exp_f1_IS_fixed_deriv_lwr), show.legend = FALSE, size= 0.5)+
  geom_ribbon(data =df_full_can %>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
                filter(sample_date >= min_date & sample_date <= max_date),
              aes(x = sample_date, ymax = ave_exp_f1_IS_fixed_deriv_upr, ymin =ave_exp_f1_IS_fixed_deriv_lwr), alpha = 0.1, size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)'")), breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 18))

ggsave(filename = "./Xiaotian et al/phac/formanuscript_nowcast_nowave/mubar_prime_ourmodel.jpeg",
       plot = gg2, 
       device = "jpeg",
       width = 5, 
       height = 3,
       dpi = 300)

gg3 <- df_full_can %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, post_prob_ave_exp_f1_IS_fixed_deriv, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_line(data =df_full_can%>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
              filter(sample_date >= min_date & sample_date <= max_date),
            aes(x = sample_date, y=post_prob_ave_exp_f1_IS_fixed_deriv), size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = "Probability", breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 14))

ggsave(filename = "./Xiaotian et al/phac/formanuscript_nowcast_nowave/prob_ourmodel.jpeg",
       plot = gg3, 
       device = "jpeg",
       width = 5, 
       height = 3,
       dpi = 300)


gg4 <- df_full_can_xiaotian %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, p1_exp_increase, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_line(data =df_full_can_xiaotian %>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
              filter(sample_date >= min_date & sample_date <= max_date),
            aes(x = sample_date, y=p1_exp_increase), size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = "Probability", breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 14))

ggsave(filename = "./Xiaotian et al/phac/formanuscript_nowcast_nowave/prob_comparison.jpeg",
       plot = gg4, 
       device = "jpeg",
       width = 5, 
       height = 3,
       dpi = 300)


gg5 <- df_full_can_xiaotian %>%
  filter(model != "fullmodel") %>% 
  group_by(model) %>% 
  filter(sample_date == max(sample_date)) %>% 
  mutate(model = factor(model,levels = paste0("model",1:103))) %>% 
  ggplot(aes(sample_date, ave_exp_f1_IS_fixed, col = model))+ 
  geom_point(size = 0.5, show.legend = FALSE)+ 
  geom_errorbar(aes(ymax = ave_exp_f1_IS_fixed_upr, ymin = ave_exp_f1_IS_fixed_lwr), show.legend = FALSE, size= 0.5)+
  geom_ribbon(data =df_full_can_xiaotian %>% filter(model == "fullmodel") %>% dplyr::select(-"model") %>% 
                filter(sample_date >= min_date & sample_date <= max_date),
              aes(x = sample_date, ymax = ave_exp_f1_IS_fixed_upr, ymin =ave_exp_f1_IS_fixed_lwr), alpha = 0.1, size = 0.2,inherit.aes = FALSE, col = "black", fill = "black")+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)")), breaks = scales::pretty_breaks(n=5),limits = c(0,200))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18))

ggsave(filename = "./Xiaotian et al/phac/formanuscript_nowcast_nowave/mubar_comparison.jpeg",
       plot = gg5, 
       device = "jpeg",
       width = 5, 
       height = 3,
       dpi = 300)


library(cowplot)
fest1 = cowplot::plot_grid(add_sub(gg1,"a) Proposed model, station average signal",size = 12),
                           add_sub(gg3,"c) Proposed model, probability of increase",size = 12),
                           add_sub(gg2,"e) Proposed model, derivative of station average signal",size = 12),
                           add_sub(gg5,"b) Comparison model, station average signal",size = 12),
                           add_sub(gg4,"d) Comparison model, probability of increase",size = 12),
                           ncol=2, align="v", byrow = FALSE) 

ggsave(filename = paste0("./Xiaotian et al/phac/formanuscript_nowcast_nowave/comparebothmodels.jpeg"),
       plot = grid.arrange(fest1), 
       device = "jpeg",
       width = 8, 
       height = 8,
       dpi = 300)




##########

### First get the plot the full canada-averaged trends both models

### Common time trend

gg1 <- df_full_can_firstwave %>% 
  filter(model == "fullmodel") %>% 
  ggplot(aes(x = sample_date, ave_exp_f1_IS_fixed)) +
  geom_line(size=0.2) + 
  geom_ribbon(aes(ymax = ave_exp_f1_IS_fixed_upr, ymin = ave_exp_f1_IS_fixed_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = "Averaged Concentrations", breaks = scales::pretty_breaks(n=5), limits = c(0,320))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



gg2 <- df_full_can %>% 
  filter(model == "fullmodel") %>% 
  ggplot(aes(x = sample_date,ave_exp_f1_IS_fixed_deriv)) + 
  geom_line(size=0.2) + 
  geom_ribbon(aes(ymax = ave_exp_f1_IS_fixed_deriv_upr, ymin = ave_exp_f1_IS_fixed_deriv_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = "Change in Averaged Concentrations", breaks = scales::pretty_breaks(n=5))+ 
  geom_hline(yintercept=0, lty="dashed")+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

fest3 = cowplot::plot_grid(gg1,gg2,ncol=1, align="v") 


gg1 <- df_full_can_xiaotian %>% 
  filter(model == "fullmodel") %>% 
  ggplot(aes(x = sample_date, ave_exp_corr_f1_IS_fixed)) +
  geom_line(size=0.2) + 
  geom_ribbon(aes(ymax = ave_exp_corr_f1_IS_fixed_upr, ymin = ave_exp_corr_f1_IS_fixed_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = "Averaged Concentrations", breaks = scales::pretty_breaks(n=5), limits = c(0,320))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



gg2 <- df_full_can_xiaotian %>% 
  filter(model == "fullmodel") %>% 
  ggplot(aes(x = sample_date,p1_exp_corr_increase)) + 
  geom_line(size=0.2) + 
  geom_point(size=0.2) + 
  theme_bw()+
  scale_y_continuous(name = "Probability Increase\n of Averaged Concentrations", breaks = scales::pretty_breaks(n=5))+ 
  geom_hline(yintercept=0.5, lty="dashed")+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "")+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

fest4 = cowplot::plot_grid(gg1,gg2,ncol=1, align="v") 

gg <- grid.arrange(fest3, fest4, nrow = 1)


ggsave(filename = "./Xiaotian et al/phac/formanuscript_phac/comparison_averaged_fullmodels.jpeg",
       plot = gg, 
       device = "jpeg",
       width = 10, 
       height = 5,
       dpi = 300)
