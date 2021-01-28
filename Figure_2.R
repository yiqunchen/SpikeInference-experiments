library(tidyverse)
library(latex2exp)
library(ggpubr)
library(SpikeInference)
plot_output_dir <- "~/Desktop/plot_output/"

construct_v_norm <- function(thj, window_size, gam) {
  tL <- max(1, thj-window_size+1)
  tR <- min(n, thj+window_size)
  v_norm <- ((gam ^ 2 * (gam ^ 2 - 1)) / (gam^2 - gam ^ (2 * (tL - thj))) + 
               (gam ^ 2 - 1) / (gam ^ (2 * (tR - thj)) - 1))
  return(v_norm)
}

construct_yphi <- function(y, thj, window_size, gam, phi) {
  yphi <- y
  tL <- max(1, thj-window_size+1)
  tR <- min(n, thj+window_size)
  v_norm <- construct_v_norm(thj, window_size, gam)
  vTy <- sum(SpikeInference::construct_v(length(y), thj, window_size, gam) * y)
  yphi[tL:thj] <- y[tL:thj] + 
    (gam * (gam ^ 2 -1) * gam ^ (- (thj - tL):0)) / (gam ^2 - gam ^ (2 * (tL - thj))) * 
    (vTy - phi) / v_norm
  yphi[(thj + 1):tR] <- y[(thj + 1):tR] - 
    (gam ^ 2 - 1) * (gam ^ (0:(tR - (thj + 1)))) / (gam ^ (2 * (tR - thj)) - 1) * 
    (vTy - phi) / v_norm
  return(yphi)
}

gam <- 0.98
n <- 80
cts <- rep(3 * gam ^(0:(n - 1)), 2)
thj <- n
h <- 5
tL <- thj - h + 1
tR <- thj + h 
phis <- seq(0, 5, length.out = 10)
firing_rate <- 0.01
target_rate <- 0.01
sigma <- 0.1

curr_sim <- simulate_ar1(n = n, gam = gam, poisMean = firing_rate,
                         sd = sigma, seed = 1)

if(sigma>=0.5){
  lam_max = 5
}else{
  lam_max = 1
}

fit_spike <- estimate_spike_by_spike_number(curr_sim, 
                                            decay_rate = gam, target_firing_rate = target_rate, 
                                            lam_min = 1e-7, lam_max = lam_max, max_iters=10,
                                            tolerance=max(5,floor(n*firing_rate*0.05)))


stopifnot(all(fit_spike$spikes==curr_sim$spikes))

thj <- fit_spike$spikes
h<- 40
nu_thj <- SpikeInference::construct_v(n, thj, h, gam)
nu_T_y <- sum(nu_thj*curr_sim$fl)
y_phi_original <- construct_yphi(curr_sim$fl, thj, h, gam, phi=nu_T_y)
# check the perturbation at nu_T_y gives us the original value
stopifnot(all(fit_spike$spikes==curr_sim$spikes))

y_phi_0 <- construct_yphi(curr_sim$fl, thj, h, gam, phi=0)
y_phi_2 <- construct_yphi(curr_sim$fl, thj, h, gam, phi=2)

df_y_phi_original <- data.frame(x=c(1:length(y_phi_original)), y=y_phi_original,
                                c = fit_spike$estimated_calcium)

fit_spike_0 <- spike_estimates(y_phi_0, decay_rate = gam, 
                               tuning_parameter = fit_spike$tuning_parameter)
fit_spike_2 <- spike_estimates(y_phi_2, decay_rate = gam, 
                               tuning_parameter = fit_spike$tuning_parameter)



df_y_phi_0 <- data.frame(x=c(1:length(y_phi_original)), y=y_phi_0,
                         c = fit_spike_0$estimated_calcium)
df_y_phi_2 <- data.frame(x=c(1:length(y_phi_original)), y=y_phi_2,
                         c = fit_spike_2$estimated_calcium)

phi_1_plot <- ggplot(df_y_phi_original%>% filter(x<=80), aes(x=x,y=y)) +
  geom_point(colour='grey50',size=1)+
  geom_line(aes(x=x,y=c),colour = "blue", size=1)+
  theme_bw()+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(limits = c(-0.6,1.9), breaks = seq(-0.5,2.0,by=0.5))+
  ylab(TeX('$y\'(\\phi)$'))+
  ggtitle(TeX("Original data ($\\phi$=$\\nu^T y$ = 1.02)"))


phi_0_plot <- ggplot(df_y_phi_0 %>% filter(x<=80), aes(x=x,y=y)) +
  geom_point(colour='grey50',size=1)+
  geom_line(aes(x=x,y=c),colour = "blue", size=1)+
  theme_bw()+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(limits = c(-0.6,1.9), breaks = seq(-0.5,2.0,by=0.5))+
  ylab(TeX('$y\'(\\phi)$'))+
  ggtitle(TeX("Perturbed data ($\\phi$=0)"))


phi_2_plot <- ggplot(df_y_phi_2 %>% filter(x<=80), aes(x=x,y=y)) +
  geom_point(colour='grey50',size=1)+
  geom_line(aes(x=x,y=c),colour = "blue", size=1)+
  theme_bw()+ 
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(limits = c(-0.6,1.9), breaks = seq(-0.5,2.0,by=0.5))+
  ylab(TeX('$y\'(\\phi)$'))+
  ggtitle(TeX("Perturbed data ($\\phi$=2)"))


#### 
set_S <- spike_inference(dat = curr_sim$fl, 
                         gam,
                         tuning_parameter = fit_spike$tuning_parameter,
                         window_size = h,
                         sig = sigma,
                         return_conditioning_sets = TRUE,
                         two_sided = FALSE)


png(paste0(plot_output_dir,'figure_2_a.png'), width = 4 ,height=2,res=200,units='in')
phi_1_plot +  xlab("Timepoint")
dev.off()

png(paste0(plot_output_dir,'figure_2_b.png'), width = 4 ,height=2,res=200,units='in')
phi_0_plot +  xlab("Timepoint")
dev.off()

png(paste0(plot_output_dir,'figure_2_c.png'), width = 4 ,height=2,res=200,units='in')
phi_2_plot +  xlab("Timepoint")
dev.off()

#### plot values for phi set

set_S$conditioning_sets[[1]]

new_x_mid <- c(c(0,2.980765e-01,9.904437e-01,2)+c(2.980765e-01,9.904437e-01,2,2.5 ))/2
#c(c(-2,-0.7717681,1.453308e-01,2)+c(-0.7717681,1.453308e-01,2,2.5))/2

plot_rec <- data.frame(x=new_x_mid,
                       y=rep(1,each=4),
                       z=factor(c(2,1,1,1)),
                       w = c(diff(c(0,2.980765e-01,9.904437e-01,2)),2.5-2))


cbPalette <- c( "#56B4E9","#E69F00", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7")


phi_values <- ggplot( plot_rec, aes(xmin = x - w/2 , xmax = x + w/2 , ymin = y, ymax = y + 1)) +
  geom_rect(aes(fill = z), colour = NA)+
  scale_fill_manual(values=cbPalette)+
  theme_bw()+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  ylab("")+
  xlab(TeX('$\\phi$'))+
  ggtitle(TeX('Range of $\\phi$ where $40 \\in M(y\'(\\phi))$ and  $\\phi>0$'))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.title.x=element_text(size=12),
        title=element_text(size=15))


png(paste0(plot_output_dir,'figure_2_d.png'), width = 12.5 ,height=1.5,res=200,units='in')
print(phi_values)
dev.off()



