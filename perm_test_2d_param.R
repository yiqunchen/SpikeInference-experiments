library(SpikeInference)
library(dplyr)

helper_function_dir <- "~/Desktop/SpikeInference-experiments/"
output_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
input_dir <- "~/Desktop/SpikeInference-experiments/input_files/"

source(paste0(helper_function_dir,"spike_inf_helper.R"))

gen_bin_vector <- function(length, position){
  result_vec <- rep(0,length)
  result_vec[position] <- 1
  return(result_vec)
}

monte_carlo_distance <- function(true_spikes, estimated_pos,
                                 num_select_spikes,
                                 num_sims = 10000,
                                 seed=1234,
                                 metric = "cor"){
  
  set.seed(seed)
  
  list_of_positions <- lapply(c(1:num_sims),
                              function(x)sample(estimated_pos, 
                                                size=num_select_spikes, 
                                                replace = FALSE))  
  
  list_of_mc_spikes <- lapply(c(1:num_sims),
                              function(i)gen_bin_vector(length(true_spikes),
                                                        list_of_positions[[i]]))
  
  if (metric=='cor'){
    list_cor_dist <- lapply(c(1:num_sims),
                            function(i)corr_metric(true_spikes, list_of_mc_spikes[[i]], 4))
  }
  
  if (metric=='vp'){
    list_cor_dist <- lapply(c(1:num_sims),
                            function(i)VictorPurpuraDist_cpp(true_spikes,
                                                             list_of_mc_spikes[[i]], 10))
  }
  
  
  if (metric=='vr'){
    list_cor_dist <- lapply(c(1:num_sims),
                            function(i)vanRossumDist_cpp(true_spikes,
                                                         list_of_mc_spikes[[i]], 0.1))
  }
  
  
  return(list("metric" = unlist(list_cor_dist)))
  
}


load(paste0(input_dir,"best_parameter_Chen_data.RData"))
h_select_seq <- c(5,20,50)
metric_name <- c('cor','vp')
exp_num_list <- c(7,8)
cell_num_array <- c(11,21,13,6,9,9,37,21,5,27)
gamma_spike_finder <- c(1-(1/100/1.25),1-(1/100/1.25),1-(1/100/2),1-(1/100/1.25),
                        1-(1/100/2),1-(1/100/1.25),1-(1/100/0.7), 1-(1/100/2),
                        1-(1/100/2),1-(1/100/0.7))
sampling_freq <- 100


for (h_select in h_select_seq){
  for (metric_name in metric_list){
  load(paste0(input_dir,"intercept_2d_tuning_param_h_",
              h_select, "_trunc_p0_mu_test_FALSE_2_sided_FALSE_large_var.RData"))
  
  df_analysis_real_data <- dplyr::bind_rows(purrr::compact(p_val_distribution_list))
  
  load(file = paste0(output_dir,"spike_2d_param_refine_results_h_",h_select,".RData"))
  df_distance <- dplyr::bind_rows(output_df_list)
  

  
  p_val_list <- vector("list",length=sum(cell_num_array[exp_num_list]))
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
      
      
      fit_calcium <- spike_estimates(calcium_1[1:T_end,curr_cell_num]+beta_0_star,
                                     gam_star, curr_lambda)$estimated_calcium
      
      test_positions <- spike_inference(calcium_1[1:T_end,curr_cell_num]+beta_0_star,
                                        gam_star, curr_lambda, window_size=h_select)
      
      current_cell <- data.frame(calcium = calcium_1[1:T_end,curr_cell_num]+beta_0_star, 
                                 spike = spike_1[1:T_end,curr_cell_num],
                                 fit_calcium=fit_calcium)
      
      current_cell <- current_cell %>% 
        mutate(cam_time = 1/100*c(1:nrow(current_cell))) 
      
      current_cell$ell_0_spikes <- 0
      current_cell[test_positions$spikes,'ell_0_spikes'] <- 1
      
      ### use permutation
      bad_name_curr_exp_num <- curr_exp_num
      bad_name_curr_cell_num <- curr_cell_num
      num_select_spikes <- df_distance %>% 
        filter(curr_exp_num==bad_name_curr_exp_num,
               curr_cell_num==bad_name_curr_cell_num) %>%
        pull(num_selective) %>% unique()
      
      
      df_metric <- df_distance %>% 
        filter(curr_exp_num==bad_name_curr_exp_num,
               curr_cell_num==bad_name_curr_cell_num,
               estimator_type=="selective") 
      
      selective_cor <- df_metric %>%
        filter(metric_type=="Cor") %>% pull(distance)
      selective_vr <- df_metric %>%
        filter(metric_type=="VR") %>% pull(distance)
      selective_vp <- df_metric %>%
        filter(metric_type=="VP") %>% pull(distance)
      
      estimated_pos <- fit_spike$spikes
      
      true_spikes <- spike_1[1:T_end,curr_cell_num] 
      cat("num_select_spikes",num_select_spikes,"\n")
      cat("estimated_pos",length(estimated_pos),"\n")
      mc_result <- monte_carlo_distance(true_spikes, estimated_pos,
                                        num_select_spikes,
                                        num_sims = 2000,
                                        seed=1234,
                                        metric = metric_name)
      
      p_metric <- mean(mc_result[['metric']] <= selective_cor)
      
      metric_ci <- quantile(mc_result[['metric']], probs = c(0.025,0.975))

      result_df <-data.frame(curr_exp_num=curr_exp_num,
                             curr_cell_num=curr_cell_num,
                             metric_lcb = metric_ci[1],
                             metric_ucb = metric_ci[2],
                             metric_name = metric_name,
                             p_metric = p_metric)
      
      p_val_list[[i]] <- list(result_df,(mc_result[['metric']]))
      i <- i+1
   
    }
  }
  
  save(p_val_list,
       file = paste0(output_dir,"new_exp_",exp_num_list,"_permutation_test_2d_param_h_",h_select,"_metric_name_",metric_name,".RData"))
  
  }
  
}
