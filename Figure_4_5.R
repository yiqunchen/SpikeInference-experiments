library(tidyverse)
library(latex2exp)
library(ggpubr)
library(RColorBrewer)
library(shades)

input_dir <- "~/Desktop/input_files/"
plot_output_dir <- "~/Desktop/plot_output/"


sim_files_h_1 <- list.files(path = input_dir,
                            pattern = glob2rx("power_mean_0.01_*h_1.*"),full.names = T)

sim_files_h_2 <- list.files(path = input_dir,
                            pattern = glob2rx("power_mean_0.01_*h_2.*"),full.names = T)

sim_files_h_10 <- list.files(path = input_dir,
                             pattern = glob2rx("power_mean_0.01_*h_10.*"),full.names = T)

sim_files_h_20 <- list.files(path = input_dir,
                             pattern = glob2rx("power_mean_0.01_*h_20.*"),full.names = T)


power_plot_01_h1 <- dplyr::bind_rows(
  lapply(sim_files_h_1,  function(x){load(x)
    dplyr::bind_rows(get("df_analysis")) }
  )) %>% mutate(h=1)

power_plot_01_h2 <- dplyr::bind_rows(
  lapply(sim_files_h_2,  function(x){load(x)
    dplyr::bind_rows(get("df_analysis")) }
  )) %>% mutate(h=2)

power_plot_01_h10 <- dplyr::bind_rows(
  lapply(sim_files_h_10,  function(x){load(x)
    dplyr::bind_rows(get("df_analysis")) }
  )) %>% mutate(h=10)

power_plot_01_h20 <- dplyr::bind_rows(
  lapply(sim_files_h_20,  function(x){load(x)
    dplyr::bind_rows(get("df_analysis")) }
  )) %>% mutate(h=20)


power_plot_01_all_h <- dplyr::bind_rows(power_plot_01_h1,power_plot_01_h2,
                                        power_plot_01_h10, power_plot_01_h20)



plot_p_val_power <- power_plot_01_all_h %>% 
  filter((pval_vec<1)) %>%
  mutate(cond_power = rejection_pow/detection_prob) %>%
  group_by(sigma_sim,h) %>%
  summarize(avg_power = mean(cond_power), 
            se_power = sd(cond_power)/sqrt(n_distinct(exp_id)),
            avg_detect_p = mean(detection_prob),
            se_detect_p = sd(detection_prob)/sqrt(n_distinct(exp_id))) %>%
  mutate(one_over_sigma = 1/sigma_sim, h = as.factor(h))


####


p_power <- plot_p_val_power %>%
  filter(one_over_sigma <= 10) %>%
  ggplot(aes(x=one_over_sigma,y=avg_power, group = h, colour = h)) + 
  geom_line(size=1)+
  geom_point(size=1)+
  ylab('Conditional power')+
  xlab(TeX("$1/\\sigma$"))+
  ggtitle('Comparison of conditional power with different SNR')+
  scale_y_continuous(breaks=seq(from=0,to=1,by=0.1))+
  labs(colour = "h")+
  theme_bw() +
  scale_color_brewer(palette="YlOrRd")+
  theme(plot.title = element_text(hjust = 0.5,size=18),
        #text = element_text(family="Times New Roman"),
        legend.position="bottom",
        axis.text = element_text(size=15),
        legend.title=element_text(size=15),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=5)))

png(paste0(plot_output_dir,'Figure_4_b','.png'), width = 6,height = 6, res=200,units='in')
p_power
dev.off()


p_detect <- plot_p_val_power %>%
  filter(one_over_sigma <= 10) %>%
  filter(h==1) %>%
  ggplot(aes(x=one_over_sigma,y=avg_detect_p),colour="black") + 
  geom_line(size=1)+
  geom_point(size=1)+
  ylab('Detection probability')+
  xlab(TeX("$1/\\sigma$"))+
  ggtitle('Detection probability')+
  #scale_y_continuous(limits=c(0.5,1),breaks=seq(from=0,to=1,by=0.1))+
  labs(colour = "Firing rate")+
  theme_bw() +
  scale_color_brewer(palette="YlOrRd")+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        #text = element_text(family="Times New Roman"),
        legend.position="bottom",
        axis.text = element_text(size=20),
        legend.title=element_text(size=20),
        legend.text = element_text(size=20,hjust = 0),
        axis.title=element_text(size=20))+
  guides(colour = guide_legend(override.aes = list(size=5)))


png(paste0(plot_output_dir,'Figure_4_a','.png'), width = 6,height = 5.5, res=200,units='in')

p_detect

dev.off()


plot_p_val_CI <- power_plot_01_all_h %>% filter(!is.na(lcb_vec) & !is.na(ucb_vec)) %>%
  mutate(cover_true_param = ((lcb_vec<=vTc)&(ucb_vec>=vTc))) %>%
  group_by(sigma_sim,h) %>%
  summarize(avg_coverage = mean(cover_true_param), 
            se_coverage = sqrt(avg_coverage*(1-avg_coverage))/sqrt(n_distinct(exp_id)),
            avg_naive_covergae = mean((naive_lcb<=vTc)&(naive_ucb>=vTc)),
            se_naive_covergae = sqrt(avg_naive_covergae*(1-avg_naive_covergae))/sqrt(n_distinct(exp_id)),
            avg_width = mean(ucb_vec-lcb_vec),
            se_width = sd(ucb_vec-lcb_vec)/sqrt(n_distinct(exp_id))) %>%
  mutate(snr = round(1/sigma_sim,1), h = as.factor(h))

### plot CI width

p_ci_coverage <- ggplot(plot_p_val_CI %>% filter(snr <= 6)
                        , aes(x=as.factor(snr), y=avg_coverage, 
                              color = h, fill = h)) + 
  geom_point(size=2, position = position_dodge(width = 0.3))+
  geom_errorbar(aes(ymin=avg_coverage-se_coverage, ymax=avg_coverage+se_coverage),
                colour="black", width=.2, position = position_dodge(width = 0.3))+
  scale_y_continuous(limits = c(0.8,1))+ #, breaks=seq(0,1,0.1))+
  theme_bw()+
  xlab(TeX("$1/\\sigma$"))+
  ylab("Coverage")+
  ggtitle(paste0('Coverage of Selective CI'))+
  labs(colour = "h")+
  theme_bw() +
  scale_color_manual(values=saturation(brewer.pal(4,"YlOrRd"), delta(0.5)))+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        #text = element_text(family="Times New Roman"),
        axis.text = element_text(size=20),
        legend.position="bottom",
        legend.title=element_text(size=20),
        legend.text = element_text(size=20,hjust = 0),
        axis.title=element_text(size=20))+
  guides(colour = guide_legend(override.aes = list(size=5)))


p_ci_width <- plot_p_val_CI %>% 
  mutate(h=as.factor(h))%>%
  filter(snr <= 6) %>%  
  mutate(sigma_sim = as.numeric(as.character(sigma_sim))) %>%
  ggplot(aes(x=1/sigma_sim,y=avg_width,group=h,
             colour=h))+
  geom_line(size=1)+
  geom_point(size=1)+
  ylab('Width of CI')+
  xlab(TeX("$1/\\sigma$"))+
  theme_bw()+
  ggtitle(paste0('Width of CIs for h=1'))+
  scale_color_manual(values=saturation(brewer.pal(4,"YlOrRd"), delta(0.5)))+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        #text = element_text(family="Times New Roman"),
        axis.text = element_text(size=20),
        legend.position="bottom",
        legend.title=element_text(size=20),
        legend.text = element_text(size=20,hjust = 0),
        axis.title=element_text(size=20))+
  guides(colour = guide_legend(override.aes = list(size=5)))



p_naive_ci_coverage <- ggplot(plot_p_val_CI %>% filter(snr <= 6)
                              , aes(x=as.factor(snr), y=avg_naive_covergae , 
                                    color = h, fill = h)) + 
  geom_point(size=2, position = position_dodge(width = 0.3))+
  geom_errorbar(aes(ymin=avg_naive_covergae-se_naive_covergae, ymax=avg_naive_covergae+se_naive_covergae),
                colour="black", width=.2, position = position_dodge(width = 0.3))+
  scale_y_continuous(limits = c(0.7,1))+ #, breaks=seq(0,1,0.1))+
  theme_bw()+
  xlab(TeX("$1/\\sigma$"))+
  ylab("Coverage")+
  ggtitle(paste0('Coverage of naive CI'))+
  labs(colour = "h")+
  theme_bw() +
  scale_color_manual(values=saturation(brewer.pal(4,"YlOrRd"), delta(0.5)))+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        #text = element_text(family="Times New Roman"),
        axis.text = element_text(size=20),
        legend.position="bottom",
        legend.title=element_text(size=20),
        legend.text = element_text(size=20,hjust = 0),
        axis.title=element_text(size=20))+
  guides(colour = guide_legend(override.aes = list(size=5)))


plot_PI_point <- power_plot_01_all_h %>% filter(!is.na(lcb_vec) & !is.na(ucb_vec)) %>%
  group_by(sigma_sim,h) %>%
  summarize(
    mid_diff_naive = mean((naive_ucb+naive_lcb)/2-vTy_vec),
    mid_diff_select =  mean((ucb_vec+lcb_vec)/2-vTy_vec ),
    se_diff_naive = sd((naive_ucb+naive_lcb)/2-vTy_vec)/sqrt(n_distinct(exp_id)),
    se_diff_select =  sd((ucb_vec+lcb_vec)/2-vTy_vec )/sqrt(n_distinct(exp_id)) ) %>%
  mutate(snr = round(1/sigma_sim,1), h = as.factor(h)) %>%
  pivot_longer(cols=c(mid_diff_select,mid_diff_naive,se_diff_select,se_diff_naive),
               names_to="type",values_to="value")


p_skew <- plot_PI_point %>% 
  filter(type%in%c("mid_diff_select"))%>%
  filter(snr <= 6) %>% 
  mutate(h=as.factor(h))%>%
  mutate(sigma_sim = as.numeric(as.character(sigma_sim))) %>%
  ggplot(aes(x=1/sigma_sim,y=value,group=interaction(h,type),
             colour=h))+
  geom_line(size=1)+
  geom_point(size=1)+
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "black", size=1)+
  ylab(TeX('Midpoint'))+
  xlab(TeX(" $1/\\sigma$"))+
  theme_bw()+
  ggtitle(TeX('Midpoint of Selective CI'))+
  scale_color_manual(values=saturation(brewer.pal(4,"YlOrRd"), delta(0.5)) )+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        axis.text = element_text(size=20),
        legend.position="bottom",
        legend.title=element_text(size=20),
        legend.text = element_text(size=20,hjust = 0),
        axis.title=element_text(size=20))+
  guides(colour = guide_legend(override.aes = list(size=5)))


naive_ci_width <-  power_plot_01_all_h %>% filter(!is.na(naive_lcb) & !is.na(naive_ucb)) %>%
  mutate(naive_width = naive_ucb-naive_lcb,
         selective_width = ucb_vec-lcb_vec)  %>%
  select(c(sigma_sim,h,exp_id,naive_width,selective_width)) %>%
  pivot_longer(c(selective_width,naive_width), names_to = "type", values_to = "width")%>%
  mutate(h=as.factor(h)) %>%
  mutate(sigma_sim = as.numeric(as.character(sigma_sim)))



p_naive_width_h_1 <- naive_ci_width%>%
  filter(h==1) %>%
  filter((1/sigma_sim) <= 6) %>% 
  ggplot()+
  geom_boxplot(aes(x=1/sigma_sim,y=width,fill=type,group=interaction(sigma_sim,type)),
               outlier.shape = NA)+
  scale_y_log10()+
  ylab('Width')+
  xlab(TeX("$1/\\sigma$"))+
  theme_bw()+
  ggtitle(paste0('Width of CIs for h=1'))+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        axis.text = element_text(size=20),
        legend.position="bottom",
        legend.title=element_text(size=20),
        legend.text = element_text(size=20,hjust = 0),
        axis.title=element_text(size=20))+
  scale_fill_manual(values = c("red","steelblue"),
                    labels=c("Wald CI","Selective CI"))+
  labs(fill="Type")


png(paste0(plot_output_dir, 'Figure_5_a','.png'),
    width = 6,height = 6, res=400,units='in')
p_ci_coverage
dev.off()

png(paste0(plot_output_dir, 'Figure_5_b','.png'),
    width = 6,height = 6, res=400,units='in')
p_naive_ci_coverage
dev.off()

png(paste0(plot_output_dir, 'Figure_5_c','.png'),
    width = 6,height = 6, res=400,units='in')
p_naive_width_h_1
dev.off()

png(paste0(plot_output_dir, 'Figure_5_d','.png'),
    width = 6,height = 6, res=400,units='in')
p_skew
dev.off()



