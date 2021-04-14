library(SpikeInference)
library(dplyr)

helper_function_dir <- "~/Desktop/SpikeInference-experiments/"
output_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
input_dir <- "~/Desktop/SpikeInference-experiments/input_files/"

source(paste0(helper_function_dir,"spike_inf_helper.R"))
load(paste0(input_dir, "best_parameter_Chen_data.RData"))
# vary h to get different result 
h_seq <- c(5,20,50)
cat("h",h,"\n")
two_sided <- FALSE
lower_trunc <- "p0" 
mu_star <- FALSE

fit_p_val_corr <- function(dat, decay_rate, h, best_lam, true_spike, h_spike=h, 
                           two_sided = two_sided, mu = 0, lower_trunc = 0){
  
  cat("current mu ", mu, "\n")
  cat("current lower_trunc ", lower_trunc, "\n")
  stopifnot(length(best_lam)==1)
  
  sigma_2_hat_jnfl <- JNFL_var_estimator(dat, decay_rate)^2
  sigma_2_hat_mad <- MAD_var_estimator(dat, decay_rate)^2  # conservative
  
  fit_calcium <- spike_estimates(dat=dat, decay_rate = decay_rate,
                                 tuning_parameter = best_lam)$estimated_calcium
  
  sigma_2_hat_calcium <- var(dat-fit_calcium) 
  
  cat("JNFL ", sigma_2_hat_jnfl, "MAD ", sigma_2_hat_mad,
      "CALCIUM ", sigma_2_hat_calcium, "OVERALL VARIANCE" ,var(dat),"\n")
  
  sigma_2_hat <- max(sigma_2_hat_jnfl,sigma_2_hat_mad,sigma_2_hat_calcium)
  
  selective_p_val <- spike_inference(dat = dat, decay_rate = decay_rate,
                                     tuning_parameter = best_lam, 
                                     window_size = h,
                                     sig2 = sigma_2_hat, 
                                     return_conditioning_sets = FALSE,
                                     two_sided = two_sided,
                                     return_ci=FALSE)
  
  n_cols <- 8
  result <- matrix(0, nrow = length(dat), ncol = n_cols)
  result[selective_p_val$change_pts,1] <- 1
  result[selective_p_val$change_pts,2] <- selective_p_val$pvals
  
  spike_loc <- selective_p_val$change_pts
  
  v_collection <- lapply(spike_loc, function(thj)
    construct_v(length(dat), thj, h, decay_rate))
  vty <- lapply(v_collection, function(v)sum(v*dat))
  vtv <- lapply(v_collection, function(v)sum(v*v))
  
  vtc_estimated <- lapply(v_collection,function(v)sum(v*fit_calcium))
  

  vtc_true <- lapply(selective_p_val$change_pts, 
                     function(x)sum(true_spike[max(1,x-h_spike+1):min(length(true_spike),
                                                                      x+h_spike)] ))
  
  vtc_true_relaxed <- lapply(selective_p_val$change_pts, 
                             function(x)sum(true_spike[max(1,min(x-10,x-h_spike+1)):min(length(true_spike),
                                                                                        x+h_spike)] ))
  
  p_val_result <- unlist(lapply(seq_along(vty), function(i)naive_z_test(vty[[i]], vtv[[i]],
                                                                        sigma2 = sigma_2_hat)))
  result[selective_p_val$change_pts,3] <- p_val_result
  result[,4] <- true_spike
  result[selective_p_val$change_pts,5] <- unlist(vtc_true)
  result[selective_p_val$change_pts,6] <- unlist(vtc_true_relaxed)
  
  result[selective_p_val$change_pts,7] <- unlist(vtc_estimated)
  result[selective_p_val$change_pts,8] <- unlist(vty)
  
  # CI
  #result[selective_p_val$change_pts,8] <- selective_p_val$LCB
  #result[selective_p_val$change_pts,9] <- selective_p_val$UCB
  
  colnames(result) <- c("estimated_spikes","selective_p_val",
                        "naive_p_val","true_spikes","vtc_true",
                        "vtc_true_relaxed", "estimated_vTc","vTy")
  
  result_df <- data.frame(result)
  
  fit_spike <- spike_estimates(dat=dat, decay_rate = decay_rate, 
                               tuning_parameter = best_lam)
  
  noise_vec <- (dat-fit_spike$estimated_calcium)
  
  return(list(result_df,noise_vec,spike_loc))
  
}

shift = function(x, lag) {
  switch(sign(lag)/2+1.5, dplyr::lead(x, abs(lag)), dplyr::lag(x, abs(lag)))
}

best_time_lag <- function(true_vTc, estimated_vTc, t_seq = seq(-20,20,1)){
  lagged_cor <- sapply(t_seq,
                       function(t)cor(shift(true_vTc,t),estimated_vTc, use="complete.obs"))
  names(lagged_cor) <- t_seq
  return(lagged_cor)
}

best_h_smooth <- function(true_spike, estimated_vTc, change_pts, h,h_adj, method="correlation"){
  stopifnot(method %in% c("correlation","binary_classification"))
  
  true_vTc <- rep(0, length(true_spike))
  true_spike[is.na(true_spike)] <- 0
  true_vTc_value <- sapply(change_pts, 
                           function(x)sum(true_spike[max(1,x-h-h_adj+1):min(length(true_spike),
                                                                            x+h+floor(h_adj))] ))
  true_vTc[change_pts] <- true_vTc_value
  
  if(method=="correlation"){
    smoothed_cor <- cor(true_vTc, estimated_vTc, use="complete.obs")
  }
  
  if(method=="binary_classification"){
    confusion_mat <- table((estimated_vTc>0),(true_vTc>0))
    smoothed_cor <- -1*confusion_mat[2,1]#1-(confusion_mat[2,1]/(confusion_mat[2,1]+confusion_mat[2,2]))
  }
  
  return(smoothed_cor)
}

adj_ground_truth <- function(true_spikes, estimated_vTc, change_pts,
                             t_seq=seq(-50,50,1), h, h_seq = seq(h,20,1)){
  #if(h>=max(h_seq)){
  #  h_seq = h
  #}
  
  # first adjust time
  cor_time_adj <- best_time_lag(true_spikes,estimated_vTc,t_seq)
  t_adj <- as.numeric(names(which.max(cor_time_adj)))
  true_spikes_adj <- shift(true_spikes, t_adj)
  
  cor_h_adj <- sapply(h_seq, function(h)
    best_h_smooth(true_spikes_adj, estimated_vTc, change_pts, h, h_seq, method = "correlation"))
  
  names(cor_h_adj) <- h_seq
  h_adj <- as.numeric(names(which.max(cor_h_adj)))+h # this is the true adjustment
  
  # return the result
  true_spikes_adj[is.na(true_spikes_adj)] <- 0 
  true_vTc_result <- rep(0, length(true_spikes))
  # averaged over the mean
  true_vTc_value <- sapply(change_pts, 
                           function(x)sum(true_spikes_adj[max(1,x-h_adj+1):min(length(true_spikes_adj),
                                                                               x+h_adj)] ))
  true_vTc_result[change_pts] <- true_vTc_value
  
  result_list <- list(true_vTc_result,cor_h_adj,cor_time_adj,t_adj,h_adj)
  names(result_list) <- c("true_vTc_result","Cor H", "Cor T", "T", "h")
  return(result_list)
}


exp_num_list= c(7,8)
cell_num_array= c(11,21,13,6,9,9,37,21,20,27)

noise_list <- vector("list",length=sum(cell_num_array))
p_val_result_list <- vector("list",length=sum(cell_num_array))
p_val_distribution_list  <- vector("list",length=sum(cell_num_array))
LAMBDA_MAX <- 100

#load(paste0(input_dir,"optimal_lambda_df.RData"))
#best_lambda_df <- best_lambda_df %>% mutate(lambda=as.numeric(lambda))

load(paste0(input_dir,"c_0_h_collection.RData"))


for (h in h_seq){
  i <- 1
  for (curr_exp_num in exp_num_list){
    
    curr_cell_list <- c(1:cell_num_array[curr_exp_num])
    
    gamma_spike_finder <- c(1-(1/100/1.25),1-(1/100/1.25),1-(1/100/2),1-(1/100/1.25),
                            1-(1/100/2),1-(1/100/1.25),1-(1/100/0.7), 1-(1/100/2),
                            1-(1/100/2),1-(1/100/0.7))
    
    h_spike_finder <- rep(1,10) 
    
    sampling_freq <- 100
    gam_star <- gamma_spike_finder[curr_exp_num]
    h_star <- h_spike_finder[curr_exp_num]
    
    calcium_1 <- read.csv(paste0(input_dir,curr_exp_num,".train.calcium.csv"))
    spike_1 <- read.csv(paste0(input_dir,curr_exp_num,".train.spikes.csv"))
    
    for (curr_cell_num in curr_cell_list){
      
      cat("curr_exp_num",curr_exp_num,
          "curr_cell_num",curr_cell_num,
          "h",h,"\n")
      
      T_half = floor(sum(!is.na(calcium_1[,curr_cell_num]))/2)
      T_end = sum(!is.na(calcium_1[,curr_cell_num]))
      
      #curr_lambda <- best_lambda_df %>% 
      #  filter(exp_num==curr_exp_num, cell_num==curr_cell_num) %>%
     #   pull(lambda)
      
      #curr_lambda <- 0.5*as.numeric(curr_lambda)
     

     # if (curr_exp_num == 8){
     #    if(curr_cell_num %in% c(2,9,19,20,21)){
     #       curr_lambda <- curr_lambda*2 
     #     }
     #   } 
     # if(curr_lambda>=LAMBDA_MAX){
     #   curr_lambda<-LAMBDA_MAX
     # }
      curr_lambda <-  result_df_exp_7_8 %>% 
      dplyr::filter(exp_num==curr_exp_num, 
                    cell_num==curr_cell_num) %>%
      pull(lambda_star)
    curr_lambda <- as.numeric(curr_lambda)

    
    beta_0_star <-result_df_exp_7_8 %>% 
      dplyr::filter(exp_num==curr_exp_num, 
                    cell_num==curr_cell_num) %>%
      pull(beta_star)
    

      curr_h <- h
      mu_star_value <- unique(result_c_0_h %>% 
                          filter(exp_id==curr_exp_num, 
                            cell_id==curr_cell_num,h==curr_h) %>%
                          pull(p_75_increase))
      

      if(lower_trunc=="None"){
        lower_trunc_value = -Inf
      }else if(lower_trunc=="p0"){
        lower_trunc_value = 0
      }else if(lower_trunc=="p75"){
        lower_trunc_value = mu_star
      }

      if(mu_star){
        mu_star_value = mu_star_value
      }else{
        mu_star_value = 0
      }

      
      #cost_beta_0_vec <- rep(NA, times = length(beta_0_range))
      
      #for (j in seq_along(beta_0_range)){
      #  beta_0 <- beta_0_range[j]
      #  current_fit <- spike_estimates(calcium_1[1:T_end,curr_cell_num]+beta_0,
      #                 gam_star, curr_lambda)
      #  current_fit_cost <- current_fit$cost[length(current_fit$cost)]
      #  cost_beta_0_vec[j] <- current_fit_cost
      #}
      
      #beta_0_star <- beta_0_range[which.min(cost_beta_0_vec)]
      
      cat("curr_exp_num",curr_exp_num,
          "curr_cell_num",curr_cell_num,
          "beta_0_star",beta_0_star,"\n")
      
      # adjusted example 
      selective_example <- fit_p_val_corr(calcium_1[1:T_end,curr_cell_num]+beta_0_star, 
                                          gam_star, h, curr_lambda,
                                          spike_1[1:T_end,curr_cell_num],
                                          two_sided = two_sided,
                                          mu = mu_star_value,
                                          lower_trunc = lower_trunc_value)
      
      p_val_result_df <- selective_example[[1]]
      p_val_change_pts <- selective_example[[3]]
      
      p_adj_result <- adj_ground_truth(p_val_result_df$true_spikes, 
                                       p_val_result_df$estimated_vTc, 
                                       p_val_change_pts,
                                       t_seq=0,#seq(-1,1,1), 
                                       h, h_seq = h_star)
      
      
      p_val_distribution <- selective_example[[1]] %>%
        mutate(vtc_true_adj = p_adj_result[[1]])%>% 
        mutate(time = 1/sampling_freq*c(1:nrow(selective_example[[1]]))) %>%
        filter(estimated_spikes==1) %>%
        mutate(exp_num = curr_exp_num, 
               cell_num = curr_cell_num,
               beta_0_star = beta_0_star,
               h=h)
      
      if(two_sided){
        
        confusion_mat <- selective_example[[1]] %>%
          mutate(vtc_true_adj = p_adj_result[[1]]) %>%
          filter(estimated_spikes==1) %>%
          select(-c(naive_p_val)) %>%
          summarise(
            num_spikes = n(),
            fp = sum(abs(vtc_true_adj)<=1e-6 & (selective_p_val<=0.05)),
            fn = sum(abs(vtc_true_adj)>1e-6 & (selective_p_val>0.05)),
            tp = sum(abs(vtc_true_adj)>1e-6 & (selective_p_val<=0.05)),
            tn = sum(abs(vtc_true_adj)<=1e-6 & (selective_p_val>0.05)),
            h_adjust = p_adj_result[[5]],
            t_adjust = p_adj_result[[4]]
          ) %>%
          mutate(exp_num = curr_exp_num, cell_num = curr_cell_num, h=h)
        
      }else{
        confusion_mat <- selective_example[[1]] %>%
          mutate(vtc_true_adj = p_adj_result[[1]]) %>%
          filter(estimated_spikes==1) %>%
          select(-c(naive_p_val)) %>%
          summarise(
            num_spikes = n(),
            fp = sum(abs(vtc_true_adj)<=1e-6 & (selective_p_val<=0.05)),
            fn = sum((vtc_true_adj)>1e-6 & (selective_p_val>0.05)),
            tp = sum((vtc_true_adj)>1e-6 & (selective_p_val<=0.05)),
            tn = sum(abs(vtc_true_adj)<=1e-6 & (selective_p_val>0.05)),
            h_adjust = p_adj_result[[5]],
            t_adjust = p_adj_result[[4]]
          ) %>%
          mutate(exp_num = curr_exp_num, cell_num = curr_cell_num, h=h)
      }
      
      p_val_distribution_list[[i]] <- p_val_distribution
      p_val_result_list[[i]] <- confusion_mat
      i <- i + 1
      
    }
  }
  #### fix time lag from -10 to 10
  #### 
  output_name = paste0(output_dir,"intercept_2d_tuning_param_h_",h,"_trunc_",lower_trunc,
    "_mu_test_",mu_star,"_2_sided_",two_sided,"_large_var",".RData")

  save(p_val_result_list,p_val_distribution_list,
         file=output_name)
  
}



