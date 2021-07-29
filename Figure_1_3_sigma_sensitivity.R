###### no cp at all ######
set.seed(1995)
library(SpikeInference)
library(tidyverse)
library(latex2exp)
library(ggpubr)


helper_function_dir <- "~/Desktop/dissertation/ell_0_project/SpikeInference-experiments/"
output_dir <- "~/Desktop/dissertation/ell_0_project/SpikeInference-experiments/input_files/"
input_dir <- "~/Desktop/dissertation/ell_0_project/SpikeInference-experiments/input_files/"
plot_output_dir <- "~/Desktop/dissertation/ell_0_project/revision-files/"

# load data and plot

sim_files <- list.files(path = input_dir,
                        pattern = "estimated_sigma_Type_1_",full.names = T)
type_1_plot <- dplyr::bind_rows(
  lapply(sim_files,  function(x){load(x)
    dplyr::bind_rows(get("df_plot")) }
  ))

library(RColorBrewer)
library(shades)
# yellow orange red palette
##### plot for one h

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

png(paste0(plot_output_dir,'supp_figure_3_b','.png'), width = 6,height = 6, res=200,units='in')
print(p_selective)
dev.off()




