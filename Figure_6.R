library(SpikeInference)
library(latex2exp)
library(dplyr)
library(ggplot2)
library(ggpubr)

helper_function_dir <- "~/Desktop/SpikeInference-experiments/"
output_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
input_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
plot_output_dir <- "~/Desktop/SpikeInference-experiments/plot_output/"

sim_files <- list.files(path = input_dir, pattern = "trunc_p0_mu_test_FALSE",full.names = T)

df_analysis_real_data <- dplyr::bind_rows(lapply(sim_files,
                                                 function(x){
                                                   load(x)
                                                   dplyr::bind_rows(get("p_val_distribution_list"))
                                                 }))

df_analysis_real_data <- df_analysis_real_data %>%
  distinct()

curr_exp_num <- 7
curr_cell_num <- 29
calcium_1 <- read.csv(paste0(input_dir,curr_exp_num,".train.calcium.csv"))
spike_1 <- read.csv(paste0(input_dir,curr_exp_num,".train.spikes.csv"))

T_end = sum(!is.na(calcium_1[,curr_cell_num]))
gamma_spike_finder <- c(1-(1/100/1.25),1-(1/100/1.25),1-(1/100/2),1-(1/100/1.25),
                        1-(1/100/2),1-(1/100/1.25),1-(1/100/0.7), 1-(1/100/2),
                        1-(1/100/2),1-(1/100/0.7))


load(input_dir,"best_parameter_Chen_data.RData")

sampling_freq <- 100
gam_star <- gamma_spike_finder[curr_exp_num]

curr_lambda <-  result_df_exp_7_8 %>% 
  dplyr::filter(exp_num==curr_exp_num, 
                cell_num==curr_cell_num) %>%
  pull(lambda_star)


beta_0_star <- result_df_exp_7_8 %>% 
  filter(exp_num==curr_exp_num, 
         cell_num==curr_cell_num)  %>% 
  pull(beta_star)

fit_spike <- spike_estimates(calcium_1[1:T_end,curr_cell_num]+beta_0_star,
                             gam_star, curr_lambda)

fit_calcium <- spike_estimates(calcium_1[1:T_end,curr_cell_num]+beta_0_star,
                               gam_star, curr_lambda)$estimated_calcium

easier_view <- df_analysis_real_data %>% 
  filter(exp_num==curr_exp_num, cell_num==curr_cell_num, h==20)  %>% 
  mutate(selective_rej=selective_p_val<=0.05) %>%
  select(c(time,selective_rej,selective_p_val,vtc_true,vtc_true_relaxed, vTy, vtc_true_adj)) %>%
  mutate(FP=as.numeric((selective_rej)&(vtc_true_adj==0)) ) 

current_cell <- data.frame(calcium = calcium_1[1:T_end,curr_cell_num]+beta_0_star, 
                           spike = spike_1[1:T_end,curr_cell_num],
                           fit_calcium=fit_calcium)

current_cell <- current_cell %>% 
  mutate(cam_time = 1/100*c(1:nrow(current_cell))) 

current_cell$no_prune_time <- 0
current_cell[fit_spike$spikes,'no_prune_time'] <- 1

current_cell_joined <- current_cell %>%
  full_join(easier_view, by=c("cam_time"="time")) %>%
  full_join(easier_view_small_h, by=c("cam_time"="time")) 

consistent_color <- scales::hue_pal()(8)

df_return <- current_cell_joined

p_return_full <- df_return %>%
  ggplot(aes(x=cam_time,y=calcium))+
  geom_point(color = "grey", size = 0.6, alpha = 0.8)+
  #geom_line(color="blue",aes(x=cam_time,y=fit_calcium))+
  xlab("Time (seconds)") +
  ylab("Fluorescence")+
  coord_cartesian(ylim=c(-4, max(df_return %>%select(calcium) * 1.001)))+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_segment(data =  df_return%>%
                 filter(no_prune_time==1)%>%
                 select(cam_time), 
               aes(x = cam_time, xend = cam_time,
                   y = -0.5, yend = -1.5, color = "orange"),
               alpha = 0.75)+
  geom_segment(data =  df_return %>%
                 filter((no_prune_time==1)&!is.na(selective_p_val)&(selective_p_val<=0.05))%>%
                 select(cam_time), 
               aes(x = cam_time, xend = cam_time,
                   y = -2.0, yend = -3.0,color = "royalblue"),
               alpha = 0.75)+
  geom_segment(data = df_return%>%
                 filter(spike>0)%>%
                 select(cam_time), 
               aes(x = cam_time, xend = cam_time,
                   y = -3.5, yend = -4.5, color = "black"),
               alpha = 0.75) +
  scale_colour_manual(name = unname(TeX(c(""))),
                      values =  c("orange","royalblue","black"),
                      breaks = c("orange", "royalblue","black"),
                      labels = unname(TeX(c("\u2113_0 solution",
                                            "\u2113_0 solution, $p<0.05$, h=20",
                                            "True spikes"))))+
  guides(colour = guide_legend(override.aes = list(size = 15)))+
  theme(legend.position="bottom",
        legend.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text = element_text(size=15))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))


p_return_zoom_1 <- df_return %>%
  filter(cam_time>=50,cam_time<=100) %>%
  ggplot(aes(x=cam_time,y=calcium))+
  geom_point(color = "grey", size = 0.6, alpha = 0.8)+
  xlab("Time (seconds)") +
  ylab("Fluorescence")+
  coord_cartesian(ylim=c(-4, max(df_return %>%
                                   filter(cam_time>=50,cam_time<=100)%>%
                                   pull(calcium)) * 1.001))+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  geom_segment(data =  df_return%>%
                 filter(no_prune_time==1)%>%
                 filter(cam_time>=50,cam_time<=100) %>%
                 select(cam_time), 
               aes(x = cam_time, xend = cam_time,
                   y = -0.5, yend = -1.5, color = "orange"),
               alpha = 0.75)+
  geom_segment(data =  df_return%>%
                 filter(cam_time>=50,cam_time<=100) %>%
                 filter((no_prune_time==1)&!is.na(selective_p_val)&(selective_p_val<=0.05))%>%
                 select(cam_time), 
               aes(x = cam_time, xend = cam_time,
                   y = -2.0, yend = -3.0,color = "royalblue"),
               alpha = 0.75)+
  geom_segment(data = df_return%>%
                 filter(cam_time>=50,cam_time<=100) %>%
                 filter(spike>0)%>%
                 select(cam_time), 
               aes(x = cam_time, xend = cam_time,
                   y = -3.5, yend = -4.5, color = "black"),
               alpha = 0.75) +
  scale_colour_manual(name = unname(TeX(c(""))),
                      values =  c("orange","royalblue","black"),
                      breaks = c("orange", "royalblue","black"),
                      labels = unname(TeX(c("\u2113_0 solution",
                                            "\u2113_0 solution, $p<0.05$, h=20",
                                            "True spikes"))))+
  guides(colour = guide_legend(override.aes = list(size = 15)))+
  theme(legend.position="bottom",
        legend.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text = element_text(size=15))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 15))


p_return_zoom_2 <- df_return %>%
  filter(cam_time>=75,cam_time<=90) %>%
  ggplot(aes(x=cam_time,y=calcium))+
  geom_point(color = "grey", size = 0.6, alpha = 0.8)+
  xlab("Time (seconds)") +
  ylab("Fluorescence")+
  coord_cartesian(ylim=c(-4, max(df_return%>%
                                   filter(cam_time>=75,cam_time<=90) %>%select(calcium) * 1.001)))+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_segment(data =  df_return%>%
                 filter(no_prune_time==1)%>%
                 filter(cam_time>=75,cam_time<=90) %>%
                 select(cam_time), 
               aes(x = cam_time, xend = cam_time,
                   y = -0.5, yend = -1.5, color = "orange"),
               alpha = 0.75)+
  geom_segment(data =  df_return%>%
                 filter(cam_time>=75,cam_time<=90) %>%
                 filter((no_prune_time==1)&!is.na(selective_p_val)&(selective_p_val<=0.05))%>%
                 select(cam_time), 
               aes(x = cam_time, xend = cam_time,
                   y = -2.0, yend = -3.0,color = "royalblue"),
               alpha = 0.75)+
  geom_segment(data = df_return%>%
                 filter(cam_time>=75,cam_time<=90) %>%
                 filter(spike>0)%>%
                 select(cam_time), 
               aes(x = cam_time, xend = cam_time,
                   y = -3.5, yend = -4.5, color = "black"),
               alpha = 0.75) +
  scale_colour_manual(name = unname(TeX(c(""))),
                      values =  c("orange","royalblue","black"),
                      breaks = c("orange", "royalblue","black"),
                      labels = unname(TeX(c("\u2113_0 solution",
                                            "\u2113_0 solution, $p<0.05$, h=20",
                                            "True spikes"))))+
  guides(colour = guide_legend(override.aes = list(size = 15)))+
  theme(legend.position="bottom",
        legend.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text = element_text(size=15))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 15))


png(paste0(plot_output_dir, 'Figure_6_a.png'),
    width = 15,height = 4, res=300,units='in')
p_return_full
dev.off()

png(paste0(plot_output_dir, 'Figure_6_b.png'),
    width = 15,height = 4, res=300,units='in')
p_return_zoom_1
dev.off()

png(paste0(plot_output_dir, 'Figure_6_c.png'),
    width = 15,height = 4, res=300,units='in')
p_return_zoom_2
dev.off()




