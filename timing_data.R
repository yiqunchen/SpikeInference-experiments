# timing experiment for T=100,000, # of spikes ~1,000
library(SpikeInference)

output_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
input_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
plot_output_dir <- "~/Desktop/SpikeInference-experiments/plot_output/"


set.seed(1995)
sim_times <- 50
spike_nums <- rep(NA, times=sim_times)
estimated_spike_nums <- rep(NA, times=sim_times)
elasped_time <- rep(NA, times=sim_times)
fit_output <- vector('list',length = sim_times)

h_seq <- c(2,5,10,20,50,100)
for (h in h_seq){
  for (i in c(1:sim_times)){
    curr_sim <- simulate_ar1(n = 10000, gam = 0.95, poisMean = 0.01, sd = 0.15, seed = i)
    spike_nums[i] <- length(curr_sim$spikes)
    
    fit_spike <- spike_estimates(dat = curr_sim$fl, decay_rate = 0.95, tuning_parameter = 0.3, 
                                 functional_pruning_out = FALSE)
    
    estimated_spike_nums[i] <- length(fit_spike$change_pts)
    
    curr_time <- system.time(inf_spike <- spike_inference(dat = curr_sim$fl, decay_rate = 0.95, 
                                         tuning_parameter = 0.3, 
                                         window_size = h, 
                                         sig2 = 0.15*0.15, 
                                         return_conditioning_sets = FALSE))
    fit_output[[i]] <- inf_spike
    elasped_time[i] <- curr_time[3]

  }
  save(elasped_time, file = paste0(output_dir,'timing_',sim_times,'_reps_no_check_h_',h,'.RData'))
  
}


