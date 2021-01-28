library(SpikeInference)
library(dplyr)

helper_function_dir <- "~/Desktop/SpikeInference-experiments/"
output_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
input_dir <- "~/Desktop/SpikeInference-experiments/input_files/"

source(paste0(helper_function_dir,"spike_inf_helper.R"))
load(paste0(input_dir,"best_parameter_Chen_data.RData"))
h_select_seq <- c(5,20,50)

summarize_distance <- function(current_cell_joined){
  
  current_cell_joined$true_spikes <- as.numeric(current_cell_joined$true_spikes>=1)
  
  cor_ell_0 <- corr_metric(current_cell_joined$true_spikes, current_cell_joined$ell_0_spikes, 4)
  cor_naive <- corr_metric(current_cell_joined$true_spikes, current_cell_joined$naive_spikes, 4)
  cor_select <- corr_metric(current_cell_joined$true_spikes, current_cell_joined$selective_spikes, 4)
  cor_naive_size <- corr_metric(current_cell_joined$true_spikes, current_cell_joined$naive_spikes_size_adjust, 4)
  
  vr_ell_0 <- vanRossumDist_cpp(current_cell_joined$true_spikes, current_cell_joined$ell_0_spikes, 0.1)
  vr_naive <- vanRossumDist_cpp(current_cell_joined$true_spikes, current_cell_joined$naive_spikes, 0.1)
  vr_select <- vanRossumDist_cpp(current_cell_joined$true_spikes, current_cell_joined$selective_spikes, 0.1)
  vr_naive_size <- vanRossumDist_cpp(current_cell_joined$true_spikes, current_cell_joined$naive_spikes_size_adjust, 0.1)
  
  vp_ell_0 <- VictorPurpuraDist_cpp(current_cell_joined$true_spikes, 
                                    current_cell_joined$ell_0_spikes, 10)
  vp_naive <- VictorPurpuraDist_cpp(current_cell_joined$true_spikes, 
                                    current_cell_joined$naive_spikes, 10)
  vp_select <- VictorPurpuraDist_cpp(current_cell_joined$true_spikes, 
                                     current_cell_joined$selective_spikes, 10)
  vp_naive_size <- VictorPurpuraDist_cpp(current_cell_joined$true_spikes,
                                         current_cell_joined$naive_spikes_size_adjust, 10)
  
  
  distance_df <- data.frame(distance = c(cor_ell_0,cor_naive,cor_select,cor_naive_size,vr_ell_0,vr_naive,vr_select,vr_naive_size,
                                         vp_ell_0,vp_naive,vp_select,vp_naive_size),
                            metric_type = rep(c("Cor", "VR","VP"),each = 4),
                            estimator_type = rep(c("ell_0","naive","selective","naive_size"),times=3),
                            num_true_spikes = sum(current_cell_joined$true_spikes),
                            num_ell_0= sum(current_cell_joined$ell_0_spikes),
                            num_selective=sum(current_cell_joined$selective_spikes),
                            num_naive = sum(current_cell_joined$naive_spikes)
  )
  return(distance_df)
  
}
exp_num_list <- c(7,8)
cell_num_array <- c(11,21,13,6,9,9,37,21,5,27)
gamma_spike_finder <- c(1-(1/100/1.25),1-(1/100/1.25),1-(1/100/2),1-(1/100/1.25),
                        1-(1/100/2),1-(1/100/1.25),1-(1/100/0.7), 1-(1/100/2),
                        1-(1/100/2),1-(1/100/0.7))
sampling_freq <- 100


for (h_select in h_select_seq){
  load(paste0(input_dir, "./intercept_2d_tuning_param_h_", 
              h_select, "_trunc_p0_mu_test_FALSE_2_sided_FALSE_large_var.RData"))
  
  df_analysis_real_data <- dplyr::bind_rows(purrr::compact(p_val_distribution_list))
  

  output_df_list <- vector("list",length=sum(cell_num_array[7:8]))
  
  i <- 1
  
  for (curr_exp_num in exp_num_list){
    calcium_1 <- read.csv(paste0(input_dir,curr_exp_num,".train.calcium.csv"))
    spike_1 <- read.csv(paste0(input_dir,curr_exp_num,".train.spikes.csv"))
    
    gam_star <- gamma_spike_finder[curr_exp_num]
    curr_cell_list <- c(1:cell_num_array[curr_exp_num])
    
    for (curr_cell_num in curr_cell_list){
      
      T_end = sum(!is.na(calcium_1[,curr_cell_num]))
      
      cat("curr_exp_num ",curr_exp_num,"curr_cell_num",curr_cell_num,"\n")
      
      curr_lambda <-  result_df_exp_7_8 %>% 
        dplyr::filter(exp_num==curr_exp_num, 
                      cell_num==curr_cell_num) %>%
        pull(lambda_star)
      curr_lambda <- as.numeric(curr_lambda)
      
      
      beta_0_star <-result_df_exp_7_8 %>% 
        dplyr::filter(exp_num==curr_exp_num, 
                      cell_num==curr_cell_num) %>%
        pull(beta_star)
      
      
      fit_spike <- spike_estimates(calcium_1[1:T_end,curr_cell_num]+beta_0_star,
                                   gam_star, curr_lambda)
      
      
      sum(spike_1[1:T_end,curr_cell_num]>=1)
      length(fit_spike$spikes)
      
      easier_view <- df_analysis_real_data %>% 
        filter(exp_num==curr_exp_num, cell_num==curr_cell_num)  %>% 
        mutate(selective_rej=selective_p_val<=0.05) %>%
        select(c(time,selective_rej,selective_p_val,
                 naive_p_val,vtc_true,vtc_true_relaxed, vTy, vtc_true_adj,h))
      
      current_cell <- data.frame(calcium = calcium_1[1:T_end,curr_cell_num]+beta_0_star, 
                                 spike = spike_1[1:T_end,curr_cell_num],
                                 fit_calcium=fit_calcium)
      
      current_cell <- current_cell %>% 
        mutate(cam_time = 1/100*c(1:nrow(current_cell))) 
      
      current_cell$no_prune_time <- 0
      current_cell[fit_spike$spikes,'no_prune_time'] <- 1
      
      current_cell_joined <- current_cell %>%
        full_join(easier_view, by=c("cam_time"="time"))
      
      current_cell_joined <- current_cell_joined %>% 
        mutate( true_spikes = spike_1[1:T_end,curr_cell_num],
                ell_0_spikes = no_prune_time,
                naive_spikes = as.numeric(!is.na(naive_p_val)&(naive_p_val<=0.05)),
                selective_spikes = as.numeric(!is.na(selective_p_val)&(selective_p_val<=0.05)))
      
      
      num_selective_spikes <- sum(current_cell_joined$selective_spikes, na.rm = TRUE)
      
      same_size_naive <- current_cell_joined %>% filter(spike==1) %>%  
        mutate(naive_ranking = rank(naive_p_val, ties.method = 'first'),
               selective_ranking = rank(selective_p_val, ties.method = 'first')) %>%
        filter(naive_ranking<=num_selective_spikes) %>%
        select(cam_time)
      
      current_cell_joined$naive_spikes_size_adjust <- 0
      current_cell_joined[current_cell_joined$cam_time%in%(same_size_naive$cam_time),"naive_spikes_size_adjust"] <- 1
      
      # size adjust
      result_df <- summarize_distance(current_cell_joined) %>%
        mutate(curr_exp_num=curr_exp_num,curr_cell_num=curr_cell_num)
      
      output_df_list[[i]] <- result_df
      i <- i+1
      
    }
  }
  
  save(output_df_list,
       file = paste0(output_dir,"spike_2d_param_refine_results_h_",h_select,".RData"))
  
}






