library(SpikeInference)
library(latex2exp)
library(dplyr)
library(ggplot2)

# change curr_h to 5 and 50 to produce figures 9 ans 10, respectively
curr_h <- 20 
set.seed(1995)

output_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
input_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
plot_output_dir <- "~/Desktop/SpikeInference-experiments/plot_output/"

load(paste0(input_dir,"new_exp_7_permutation_test_2d_param_h_",curr_h,"_metric_name_vp.RData"))
exp_7_vp_result_all_pos <- dplyr::bind_rows(lapply(p_val_list,function(x)x[[1]]))

load(paste0(input_dir,"new_exp_8_permutation_test_2d_param_h_",curr_h,"_metric_name_vp.RData"))
exp_8_vp_result_all_pos <- dplyr::bind_rows(lapply(p_val_list,function(x)x[[1]]))

vp_sim_all_pos <- rbind(exp_7_vp_result_all_pos, exp_8_vp_result_all_pos)

load(paste0(input_dir,"new_exp_7_permutation_test_2d_param_h_",curr_h,"_metric_name_cor.RData"))
exp_7_cor_result_all_pos <- dplyr::bind_rows(lapply(p_val_list,function(x)x[[1]]))

load(paste0(input_dir, "new_exp_8_permutation_test_2d_param_h_",curr_h,"_metric_name_cor.RData"))
exp_8_cor_result_all_pos <- dplyr::bind_rows(lapply(p_val_list,function(x)x[[1]]))

cor_sim_all_pos <- rbind(exp_7_cor_result_all_pos, exp_8_cor_result_all_pos) 

###### compilation of results ######
load(paste0(input_dir,"spike_2d_param_refine_results_h_",curr_h,".RData"))
df_distance <- dplyr::bind_rows(output_df_list)

#### plot correlation ####

cor_perm_plot_selective <- df_distance %>% 
  filter(metric_type=="Cor", estimator_type == "selective")  %>% 
  left_join(cor_sim_all_pos, by = c("curr_exp_num","curr_cell_num")) %>%
  left_join(df_distance %>%
              filter(metric_type=="Cor", estimator_type == "ell_0") %>%
              select(c("curr_exp_num","curr_cell_num","distance")),
            by = c("curr_exp_num","curr_cell_num"))

cor_perm_plot_selective$distance_select <- cor_perm_plot_selective$distance.x
cor_perm_plot_selective$distance_ell_0 <- cor_perm_plot_selective$distance.y
cor_perm_plot_selective$distance.y <- NULL
cor_perm_plot_selective$distance.x <- NULL

cor_perm_plot_selective$recording_num <- as.numeric(rownames(cor_perm_plot_selective))

p_cor_h <- cor_perm_plot_selective %>% 
  ggplot()+
  geom_point(aes(x=recording_num, y=distance_select, color="royalblue"),size=2)+
  geom_point(aes(x=recording_num, y=distance_ell_0, color="orange"),size=2)+
  geom_linerange(aes(x=recording_num,ymin=metric_lcb, ymax=metric_ucb),
                 position=position_dodge(.9)) +
  theme_bw() + 
  xlab("Recording number")+
  ylab("Correlation")+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=15), 
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))+
  scale_colour_manual(name = unname(TeX(c(""))),
                      labels = unname(TeX(c(paste0("\u2113_0 solution, $p<0.05, h=$",curr_h),
                                            "\u2113_0 solution"))),
                      breaks=c("royalblue","orange"),
                      values=c("royalblue","orange")) 


png(paste0(plot_output_dir,'Figure_7a.png'),
    width = 12,height = 4, res=300,units='in')
ggpubr::ggarrange(p_cor_h,
                  ncol=1,common.legend = TRUE, legend = "bottom",
                  align = "hv")
dev.off()


#### plot vp ####

vp_perm_plot_selective <- df_distance %>% 
  filter(metric_type=="VP", estimator_type == "selective")  %>% 
  left_join(vp_sim_all_pos, by = c("curr_exp_num","curr_cell_num")) %>%
  left_join(df_distance %>%
              filter(metric_type=="VP", estimator_type == "ell_0") %>%
              select(c("curr_exp_num","curr_cell_num","distance")),
            by = c("curr_exp_num","curr_cell_num"))

vp_perm_plot_selective$distance_select<- vp_perm_plot_selective$distance.x
vp_perm_plot_selective$distance_ell_0 <- vp_perm_plot_selective$distance.y
vp_perm_plot_selective$recording_num <- as.numeric(rownames(vp_perm_plot_selective))

p_vp_h <- vp_perm_plot_selective %>% 
  ggplot()+
  geom_point(aes(x=recording_num, y=distance_select, color="royalblue"),size=2)+
  geom_point(aes(x=recording_num, y=distance_ell_0, color="orange"),size=2)+
  geom_linerange(aes(x=recording_num,ymin=metric_lcb, ymax=metric_ucb),
                 position=position_dodge(.9)) +
  theme_bw() + 
  xlab("Recording number")+
  ylab("Victor-Purpura Distance")+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=15), 
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))+
  scale_colour_manual(name = unname(TeX(c(""))),
                      labels = unname(TeX(c(paste0("\u2113_0 solution, $p<0.05, h=$",curr_h),
                                            "\u2113_0 solution"))),
                      breaks=c("royalblue","orange"),
                      values=c("royalblue","orange")) 



png(paste0(plot_output_dir,'Figure_7b.png'),
    width = 12,height = 4, res=300,units='in')
ggpubr::ggarrange(p_vp_h,
                  ncol=1,common.legend = TRUE, legend = "bottom",
                  align = "hv")
dev.off()







