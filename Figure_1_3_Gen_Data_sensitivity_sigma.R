###### no cp at all ######
set.seed(1995)
library(SpikeInference)
library(tidyverse)

helper_function_dir <- "~/Desktop/SpikeInference-experiments/"
output_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
input_dir <- "~/Desktop/SpikeInference-experiments/input_files/"

# helper_function_dir <- "~/Desktop/dissertation/ell_0_project/SpikeInference-experiments/"
# output_dir <- "~/Desktop/dissertation/ell_0_project/SpikeInference-experiments/input_files/"
# input_dir <- "~/Desktop/dissertation/ell_0_project/SpikeInference-experiments/input_files/"

source(paste0(helper_function_dir,"spike_inf_helper.R"))

pval_vec <- vector("list", length = sim_times)
sim_times <- 100
gam <- 0.98

naive_pval_vec <- vector("list", length = sim_times)
selective_pval_vec <- vector("list", length = sim_times)

naive_z_test <- function(vty, vtv, sigma2, mu = 0, one.sided=TRUE){
  z_score <- vty/sqrt(sigma2*vtv)
  if(one.sided){
    p_val <- pnorm(z_score, mean = mu, lower.tail = FALSE)
  }else{
    p_val <- pnorm(abs(z_score), mean = mu, lower.tail = FALSE)+
      pnorm(-abs(z_score), mean = mu, lower.tail = TRUE)
  }
  return(max(min(1,p_val),0))
}

#LAM_MAX <- 4
#LAM_MIN <- 1e-6
target_rate <- 0.01
firing_rate <- target_rate
n_length <- 10000
LAMBDA <- 0.1
sigma <- 0.2

for (h in c(1,2,10,20)){
  for (i in c(1:sim_times)){
    cat('currently at sim ', i, "h =",h, '\n')
    
    curr_sim <- simulate_ar1(n = n_length, gam = gam, poisMean = 0.0, sd = sigma, seed = i)
    # 0.3 is tuned to have 
    
    fit_spike <- spike_estimates(dat = curr_sim$fl, decay_rate = gam, tuning_parameter = LAMBDA, 
                                 functional_pruning_out = FALSE)
    
    
    positive_selective_spike <- spike_inference(dat = curr_sim$fl, decay_rate = gam,
                                                tuning_parameter = LAMBDA, window_size = h, 
                                                sig2 = NULL, 
                                                return_conditioning_sets = FALSE,
                                                return_ci = FALSE)
    
    selective_pval_vec[[i]] <- positive_selective_spike$pvals
    
    # extract spike locations
    spike_loc <- fit_spike$change_pts
    v_collection <- lapply(spike_loc, function(thj)
      construct_v(n_length, thj, h, gam))
    vty <- lapply(v_collection, function(v)sum(v*curr_sim$fl))
    vtv <- lapply(v_collection, function(v)sum(v*v))
    
    p_val_result <- unlist(lapply(seq_along(vty), function(i)naive_z_test(vty[[i]], vtv[[i]],
                                                                          sigma2 = positive_selective_spike$sig2)))
    
    
    naive_pval_vec[[i]] <- p_val_result[unlist(vty)>=0]
    
    stopifnot(length(naive_pval_vec[[i]])==length(selective_pval_vec[[i]]))
    
  }
  
  df_plot <- data.frame(selective_p_val = unlist(selective_pval_vec), 
                        naive_p_val = unlist(naive_pval_vec))
  
  df_plot <- df_plot %>% mutate(h=h) 
  save(df_plot, file = paste0(output_dir,"estimated_sigma_Type_1_error_h_",h,'.RData'))
  
}