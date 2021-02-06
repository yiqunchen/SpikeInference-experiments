library(SpikeInference)
library(dplyr)

helper_function_dir <- "~/Desktop/SpikeInference-experiments/"
output_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
input_dir <- "~/Desktop/SpikeInference-experiments/input_files/"

source(paste0(helper_function_dir,"spike_inf_helper.R"))

exp_num_list= c(7,8)
cell_num_array= c(11,21,13,6,9,9,37,21,20,27)

cost_list <- vector("list",length=sum(cell_num_array))

i <- 1

beta_0_range <- seq(-0.8,0.15,0.05)
lambda_range <- (10^seq(-2,2,0.2))


for (curr_exp_num in exp_num_list){
  curr_cell_list <- c(1:cell_num_array[curr_exp_num])
  if (curr_exp_num == 7){
    curr_cell_list <- curr_cell_list[curr_cell_list]
  }
  if (curr_exp_num == 8){
    curr_cell_list <- curr_cell_list[curr_cell_list]
  }
  
  gamma_spike_finder <- c(1-(1/100/1.25),1-(1/100/1.25),1-(1/100/2),1-(1/100/1.25),
                          1-(1/100/2),1-(1/100/1.25),1-(1/100/0.7), 1-(1/100/2),
                          1-(1/100/2),1-(1/100/0.7))
  
  sampling_freq <- 100
  gam_star <- gamma_spike_finder[curr_exp_num]

  calcium_1 <- read.csv(paste0(input_dir,curr_exp_num,".train.calcium.csv"))
  spike_1 <- read.csv(paste0(input_dir,curr_exp_num,".train.spikes.csv"))
  
  for (curr_cell_num in curr_cell_list){
    
    cat("curr_exp_num",curr_exp_num,
        "curr_cell_num",curr_cell_num,
        "\n")
    
    T_half = floor(sum(!is.na(calcium_1[,curr_cell_num]))/2)
    T_end = sum(!is.na(calcium_1[,curr_cell_num]))
    
    cost_matrix_full <- matrix(NA, nrow = length(beta_0_range), ncol=length(lambda_range))
    cost_matrix_half <- matrix(NA, nrow = length(beta_0_range), ncol=length(lambda_range))
    
    firing_matrix_full <- matrix(NA, nrow = length(beta_0_range), ncol=length(lambda_range))
    firing_matrix_half <- matrix(NA, nrow = length(beta_0_range), ncol=length(lambda_range))
    
    for (j in seq_along(beta_0_range)){
      for (k in seq_along(lambda_range)){
        curr_beta_0 <- beta_0_range[j]
        curr_lambda <- lambda_range[k]
        current_fit_full <- spike_estimates(calcium_1[1:T_end,curr_cell_num] + curr_beta_0,
                                       gam_star, curr_lambda)
        
        current_fit_half <- spike_estimates(calcium_1[1:T_half,curr_cell_num] + curr_beta_0,
                                       gam_star, curr_lambda)
     #   cat(current_fit_full$cost[length(current_fit_full$cost)],"cost",
     # length(current_fit_full$spikes),"firing","\n") 
        cost_matrix_full[j,k] <- current_fit_full$cost[length(current_fit_full$cost)]
        cost_matrix_half[j,k] <- current_fit_half$cost[length(current_fit_half$cost)]
        
        firing_matrix_full[j,k] <- length(current_fit_full$spikes)
        firing_matrix_half[j,k] <- length(current_fit_half$spikes)
       #cat(cost_matrix_full,"\n") 
      }
    }
    
    cost_list[[i]] <- list(cost_matrix_full,
      cost_matrix_half,firing_matrix_full,
      firing_matrix_half,curr_exp_num,curr_cell_num)

    i <- i + 1
  }
  
}

output_name = paste0(output_dir,"Chen_et_al_2D_tuning_parameter.RData")

save(cost_list, file=output_name)



