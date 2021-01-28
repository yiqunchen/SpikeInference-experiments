library(ggplot2)
library(dplyr)
sim_times <- 50
output_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
input_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
plot_output_dir <- "~/Desktop/SpikeInference-experiments/plot_output/"

load(input_dir,'timing_50_reps_no_check_h_2.RData')
time_h_2 <- elasped_time
load(input_dir,'timing_50_reps_no_check_h_5.RData')
time_h_5 <- elasped_time
load(input_dir,'timing_50_reps_no_check_h_10.RData')
time_h_10 <- elasped_time
load(input_dir,'timing_50_reps_no_check_h_20.RData')
time_h_20 <- elasped_time
load(input_dir,'timing_50_reps_no_check_h_50.RData')
time_h_50 <- elasped_time
load(input_dir,'timing_50_reps_no_check_h_100.RData')
time_h_100 <- elasped_time


df_plot <- data.frame(
  time = c(time_h_2,time_h_5,time_h_10,time_h_20,time_h_50,time_h_100),
  h_size = c (rep(2, length(time_h_2)), rep(5, length(time_h_5)),
            rep(10, length(time_h_10)), rep(20, length(time_h_20)), 
                rep(50, length(time_h_50)), rep(100, length(time_h_100)) ) )

df_plot$h_size <- as.factor(df_plot$h_size)


df_plot$h_size <- as.numeric(as.character(df_plot$h_size))

p2 <- ggplot(df_plot%>%filter(h_size>=2), 
             aes(x=h_size, y=time)) + 
  geom_point() + 
  geom_smooth(method = lm, formula = (y) ~ I(x^2)+I(x), se = FALSE)+
  ylab('Time(s)')+
  ggtitle('Timing experiment for T=10,000')+
  xlab('Window size: h')+
  theme_bw() +
  scale_y_log10()+
  scale_x_log10()

png(plot_output_dir,'Figure_8.png',
    width = 7 ,height = 5,res=300,units='in')
p2+
  theme(plot.title = element_text(hjust = 0.5,size=18),
        legend.title=element_text(size=15),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=15),
        legend.position="bottom")
dev.off()





