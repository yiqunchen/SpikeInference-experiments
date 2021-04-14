###### no cp at all ######
set.seed(1995)
library(SpikeInference)
library(tidyverse)
library(latex2exp)
library(ggpubr)

helper_function_dir <- "~/Desktop/SpikeInference-experiments/"
output_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
input_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
plot_output_dir <- "~/Desktop/SpikeInference-experiments/plot_output/"

source(paste0(helper_function_dir,"spike_inf_helper.R"))

###### type I error control
gam <- 0.98
LAMBDA <- 0.09
sigma <- 0.2
n_length <- 10000
h <- 1

curr_sim <- simulate_ar1(n = n_length, gam = gam, poisMean = 0.0, sd = sigma, seed = 1)

fit_spike <- spike_estimates(dat = curr_sim$fl, decay_rate = gam, tuning_parameter = LAMBDA, 
                             functional_pruning_out = FALSE)

inference_spike_toy <- spike_inference(dat = curr_sim$fl, decay_rate = gam,
                                       tuning_parameter = LAMBDA, window_size = h, 
                                       sig = sigma*sigma, return_conditioning_sets = FALSE,
                                       return_ci = FALSE, lower_trunc = 0)

spike_loc <- inference_spike_toy$spikes
v_collection <- lapply(spike_loc, function(thj)
  construct_v(n_length, thj, h, gam))
vty <- lapply(v_collection, function(v)sum(v*curr_sim$fl))
vtv <- lapply(v_collection, function(v)sum(v*v))

p_val_result_naive_toy <- unlist(lapply(seq_along(vty), 
                                        function(i)naive_z_test(vty[[i]], vtv[[i]],sigma2 = sigma*sigma)))

# select the subset with positive test stats only


plot_toy_null <- data.frame(time = c(1:n_length),
                            fluo = fit_spike$dat, 
                            curve = fit_spike$estimated_calcium,
                            spike_loc = NA,
                            selective_pval = NA,
                            naive_pval = NA
)
plot_toy_null[spike_loc,"spike_loc"] <- 1
plot_toy_null[spike_loc,"selective_pval"] <- inference_spike_toy$pvals
plot_toy_null[spike_loc,"naive_pval"] <- p_val_result_naive_toy


p_toy <- plot_toy_null %>%
  filter(time>=3000, time<3100) %>%
  ggplot()+
  geom_point(aes(x=time,y=fluo),colour='grey60',size=1.5)+
  geom_line(aes(x=time,y=curve), colour="blue",size=1)+
  theme_bw()+
  ylab("")+
  xlab("Time")+
  ggtitle(TeX("$\u2113_0$ estimated spikes under global null"))+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        legend.position="bottom",
        axis.text = element_text(size=18),
        legend.title=element_text(size=15),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=20))

png(paste0(plot_output_dir,'figure_1_a','.png'),
    width = 6,height = 6, res=200,units='in')
p_toy
dev.off()


# load data and plot

sim_files <- list.files(path = input_dir,
                        pattern = "Type_1_error_h_",full.names = T)
type_1_plot <- dplyr::bind_rows(
  lapply(sim_files,  function(x){load(x)
    dplyr::bind_rows(get("df_plot")) }
  ))

library(RColorBrewer)
library(shades)
# yellow orange red palette
##### plot for one h

png(paste0(plot_output_dir,'figure_1_c','.png'),
    width = 6,height = 6, res=200,units='in')

p_plot <- type_1_plot %>%
  mutate(h=as.factor(h)) %>%
  filter(h==1) %>%
  arrange(selective_p_val) %>%
  mutate(theoretical = c(1:nrow(.))/nrow(.)) %>%
  ggplot() +
  geom_point(aes(y = selective_p_val, x = (theoretical)),size=0.6) +
  ylab('Selective p-value Quantiles')+
  xlab('Uniform(0,1) Quantiles')+
  ggtitle(TeX('Selective test p-values'))+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5,size=20),
        legend.position="bottom",
        legend.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.text = element_text(size=20,hjust = 0),
        axis.title=element_text(size=20))

print(p_plot)

dev.off()

#####
my_palette <- brightness(brewer.pal(4,"YlOrRd"),0.9)
my_YlOrRd <- c("#FED976","#FD8D3C","#E31A1C" ,"#B10026")

p_selective <- type_1_plot %>%
  mutate(h=as.factor(h)) %>%
  group_by(h) %>%
  mutate(theoretical = ecdf(selective_p_val)(selective_p_val)) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(y = selective_p_val, x = (theoretical), colour = h),size=0.6) +
  ylab('Selective p-value Quantiles')+
  xlab('Uniform(0,1) Quantiles')+
  ggtitle(TeX('Selective p-values under global null'))+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw() +
  scale_color_manual(values=saturation(brewer.pal(4,"YlOrRd"), delta(0.5)))+
  #scale_color_brewer(palette=)+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        #text = element_text(family = "GillSans-Light"),
        legend.position="bottom",
        legend.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.text = element_text(size=20,hjust = 0),
        axis.title=element_text(size=20))+
  guides(colour = guide_legend(override.aes = list(size=5)))

png(paste0(plot_output_dir,'figure_3_b','.png'), width = 6,height = 6, res=200,units='in')
print(p_selective)
dev.off()


p_naive <- type_1_plot %>%
  mutate(h=as.factor(h)) %>%
  group_by(h) %>%
  mutate(theoretical = ecdf(naive_p_val)(naive_p_val)) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(y = naive_p_val, x = (theoretical), colour = h),size=0.6) +
  ylab('Naive p-value Quantiles')+
  xlab('Uniform(0,1) Quantiles')+
  ggtitle(TeX('Naive p-values under global null'))+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw() +
  scale_color_manual(values=saturation(brewer.pal(4,"YlOrRd"), delta(0.5)))+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        # text = element_text(family = "GillSans-Light"),
        legend.position="bottom",
        axis.text = element_text(size=20),
        legend.title=element_text(size=20),
        legend.text = element_text(size=20,hjust = 0),
        axis.title=element_text(size=20))+
  guides(colour = guide_legend(override.aes = list(size=5)))

png(paste0(plot_output_dir, 'figure_3_a','.png'), width = 6,height = 6, res=200,units='in')
print(p_naive)
dev.off()



p_plot <- type_1_plot %>%
  mutate(h=as.factor(h)) %>%
  group_by(h) %>%
  mutate(theoretical = ecdf(naive_p_val)(naive_p_val)) %>%
  ungroup() %>% 
  filter(h==1) %>%
  ggplot() +
  geom_point(aes(y = naive_p_val, x = (theoretical)),size=0.6) +
  ylab('Naive p-value Quantiles')+
  xlab('Uniform(0,1) Quantiles')+
  ggtitle(TeX('Naive test p-values'))+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=20),
        legend.position="bottom",
        legend.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.text = element_text(size=20,hjust = 0),
        axis.title=element_text(size=20))

png(paste0(plot_output_dir,'figure_1_b','.png'),
    width = 6,height = 6, res=200,units='in')
print(p_plot)
dev.off()




