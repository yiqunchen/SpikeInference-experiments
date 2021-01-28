###### no cp at all ######
set.seed(1995)
library(SpikeInference)
library(tidyverse)
helper_function_dir <- "~/Desktop/"
output_dir <- "~/Desktop/input_files/"
source(paste0(helper_function_dir,"spike_inf_helper.R"))

# try_fit <- estimate_spike_by_spike_number(curr_sim, decay_rate = gam, target_firing_rate = 0.01, 
#     lam_min = 1e-6, lam_max = 1e-1, max_iters=50, tolerance=5)
sim_times <- 100
gam <- 0.98
sign_vec <- vector("list", length = sim_times)
pval_vec <- vector("list", length = sim_times)

sigma <- 0.2
LAM_MAX <- 4
LAM_MIN <- 1e-6
target_rate <- 0.01
firing_rate <- target_rate
n_length <- 10000
LAMBDA <- 0.1

for (h in c(1,2,10,20)){
  for (i in c(1:sim_times)){
    cat('currently at sim ', i, "h =",h, '\n')
    
    curr_sim <- simulate_ar1(n = n_length, gam = gam, poisMean = 0.0, sd = sigma, seed = i)
    # 0.3 is tuned to have 
    
    fit_spike <- spike_estimates(dat = curr_sim$fl, decay_rate = gam, tuning_parameter = LAMBDA, 
                                 functional_pruning_out = FALSE)
    
    
    inference_spike <- spike_inference(dat = curr_sim$fl, decay_rate = gam,
                                       tuning_parameter = LAMBDA, window_size = h, 
                                       sig = sigma*sigma, return_conditioning_sets = FALSE,
                                       return_ci = FALSE)
    
    spike_signs <- fit_spike$spike_sign
    spike_pvals <- inference_spike$pvals
    sign_vec[[i]] <- spike_signs
    pval_vec[[i]] <- spike_pvals
    
  }
  
  df_plot <- data.frame(p_val = unlist(pval_vec), sign_vec = unlist(sign_vec))
  
  df_plot <- df_plot %>% arrange(p_val) %>% 
    mutate(h=h, 
           theoretical=c(1:nrow(df_plot))*1/nrow(df_plot))
  
  save(df_plot, file = paste0(output_dir, "type_1_error_h_", h, '.RData'))
  
}

sim_times <- 100
gam <- 0.95
sign_vec <- vector("list", length = sim_times)
pval_vec <- vector("list", length = sim_times)

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

for (h in c(1,2,10,20)){
  for (i in c(1:sim_times)){
    cat('currently at sim ', i, "h =",h, '\n')
    
    curr_sim <- simulate_ar1(n = n_length, gam = gam, poisMean = 0.0, sd = sigma, seed = i)
    # 0.3 is tuned to have 
    fit_spike <- spike_estimates(dat = curr_sim$fl, decay_rate = gam, tuning_parameter = LAMBDA, 
                                 functional_pruning_out = FALSE)
    
    spike_loc <- fit_spike$change_pts
    v_collection <- lapply(spike_loc, function(thj)
      construct_v(n_length, thj, h, gam))
    vty <- lapply(v_collection, function(v)sum(v*curr_sim$fl))
    vtv <- lapply(v_collection, function(v)sum(v*v))
    
    p_val_result <- unlist(lapply(seq_along(vty), function(i)naive_z_test(vty[[i]], vtv[[i]],
                                                                          sigma2 = sigma*sigma)))
    
    spike_signs <- fit_spike$spike_sign
    spike_pvals <- p_val_result
    sign_vec[[i]] <- spike_signs
    pval_vec[[i]] <- spike_pvals
    
  }
  
  df_plot <- data.frame(p_val = unlist(pval_vec), sign_vec = unlist(sign_vec))
  
  df_plot <- df_plot %>% arrange(p_val) %>% 
    mutate(h=h, 
           theoretical=c(1:nrow(df_plot))*1/nrow(df_plot))
  
  save(df_plot, file = paste0(output_dir, "naive_type_1_error_h_", h, '.RData'))
  
}