library(Rcpp)


nearest_changepoint <- function(pt, cps) {
  M <- length(pt)
  out <- numeric(M)
  for (i in 1:M) {
    out[i] <- cps[which.min(abs(pt[i] - cps))]
  }
  return(out)
}

loss_function_n_spikes <- function(lam, gcamp, decay_rate, num_spikes_target){
  fit <- spike_estimates(dat = gcamp, decay_rate = decay_rate, tuning_parameter = lam, 
                         functional_pruning_out = FALSE)
  num_spiked_spike <- length(fit$spikes)
  spike_num_diff <- abs(num_spiked_spike-num_spikes_target)
  return(spike_num_diff)
}

one_d_binary_search <- function(gcamp, decay_rate, lam_min, lam_max, 
                                num_spikes_target, max_iters=50, tolerance=5){
  iter_i <- 0
  loss_i <- Inf
  verbose <- 0
  
  while(iter_i <= max_iters & loss_i > tolerance){
    
    lam_1 <- (3 * lam_min + lam_max) / 4
    lam_2 <- (3 * lam_max + lam_min) / 4
    
    loss_lam_1 = loss_function_n_spikes(lam_1, gcamp, decay_rate, num_spikes_target)
    loss_lam_2 = loss_function_n_spikes(lam_2, gcamp, decay_rate, num_spikes_target)
    
    if(verbose){
      cat('at iteration i = ', iter_i, "\n")
      cat("loss at lam 1 (" , lam_1 , ") = " , loss_lam_1, "\n")
      cat("loss at lam 2 (" , lam_2 , ") = " , loss_lam_2, "\n")
    }
    
    if (loss_lam_1 < loss_lam_2){
      if (loss_lam_1 < tolerance){
        return(lam_1)
      }
      lam_max = (lam_min + lam_max) / 2
    }else{
      if (loss_lam_2 < tolerance){
        return(lam_2)
      }
      lam_min = (lam_min + lam_max) / 2
    }
      
    iter_i = iter_i + 1
  }
  
  lam_star = (lam_min + lam_max) / 2
  
  return(lam_star)
}

estimate_spike_by_spike_number <- function(dat, decay_rate, target_firing_rate, 
                                           lam_min = 1e-6, lam_max = 1, max_iters=50, tolerance=5){
  
  fps = 1 # sampling ratio - default to one
  gcamp = dat$fl
  # transform target firing rate into a number of spikes
  n = length(gcamp)
  nbin_spike_target = floor(target_firing_rate * n / fps)
  lam_star = one_d_binary_search(gcamp, decay_rate, lam_min, lam_max,
                                  nbin_spike_target, max_iters, tolerance)
  fit = spike_estimates(gcamp, decay_rate, lam_star, 
                        functional_pruning_out = FALSE)
  return(fit)
  
}

MAD_var_estimator <- function(y, decay_rate){
  lag_1_diff <- (y[2:length(y)]-decay_rate*y[1:(length(y)-1)])/sqrt(2)
  sigma_hat <- stats::mad(lag_1_diff)
  return(sigma_hat)
}

construct_v <- function(n, thj, window_size, gam) {
  tL <- max(1, thj-window_size+1)
  tR <- min(n, thj+window_size)
  v <- 0 * numeric(n)
  thj <- max(thj,1)
  if(length(-(thj - tL):0)!=(thj-tL+1)
     ){
    stop(paste0("wrong size!",n," ", thj," ", window_size," ", gam,"\n"))
  }
  v[tL:thj] <- -gam * (gam ^ 2 - 1) / (gam ^ 2 - gam ^ (2 * (tL - thj) ) ) * gam ^ (-(thj - tL):0)
  v[(thj + 1):tR] <- (gam ^ 2 - 1) / (gam ^ (2 * (tR - thj)) - 1) * gam ^ (0:(tR - (thj + 1)))
  return(v)
}

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

fit_test_set <- function(dat, decay_rate, h, best_lam, true_spike){
  
  stopifnot(length(best_lam)==8)
  
  sigma_2_hat <- JNFL_var_estimator(dat, decay_rate)^2
  
  unconstrained_spike <- spike_estimates(dat, decay_rate, best_lam[4])
  positive_spike <- spike_estimates(dat, decay_rate, best_lam[5])
  naive_p_result <- spike_naive_test(dat, decay_rate, best_lam[3], sigma_2_hat, window_size=h)
  naive_p_result_Bon <- spike_naive_test(dat, decay_rate, best_lam[2], sigma_2_hat, window_size=h)
  naive_p_result_BH <- spike_naive_test(dat, decay_rate, best_lam[1], sigma_2_hat, window_size=h)
  
  selective_p_val <- spike_inference(dat = dat, decay_rate = decay_rate,
                                     tuning_parameter = best_lam[6], 
                                     window_size = h,
                                     sig2 = sigma_2_hat, 
                                     return_conditioning_sets = FALSE,
                                     return_ci=FALSE)

  selective_p_val_BH <- spike_inference(dat = dat, decay_rate = decay_rate,
                                     tuning_parameter = best_lam[7], 
                                     window_size = h,
                                     sig2 = sigma_2_hat, 
                                     return_conditioning_sets = FALSE,
                                     return_ci=FALSE)

  p_selective_p_val_BH <- p.adjust(selective_p_val_BH$pvals, "BH")


  selective_p_val_Bon <- spike_inference(dat = dat, decay_rate = decay_rate,
                                     tuning_parameter = best_lam[8], 
                                     window_size = h,
                                     sig2 = sigma_2_hat, 
                                     return_conditioning_sets = FALSE,
                                     return_ci=FALSE)

  p_selective_p_val_Bon <- p.adjust(selective_p_val_BH$pvals,"bonferroni")

  result_colnames <- c("normal","positive","naive","Bon","BH","selective",
  "selective_BH","selective_Bon", "true", "calcium")
  n_cols <- 10
  result <- matrix(0, nrow = length(dat), ncol = n_cols)
  colnames(result) <- result_colnames
  
  result[unconstrained_spike$change_pts,1] <- 1
  result[positive_spike$change_pts[positive_spike$spike_sign=="Positive"],2] <- 1
  result[naive_p_result$spike_loc[naive_p_result$p_val<=0.05],3] <- 1
  result[naive_p_result_Bon$spike_loc[naive_p_result_Bon$p_bonferroni<=0.05],4] <- 1
  result[naive_p_result_BH$spike_loc[naive_p_result_BH$p_BH<=0.05],5] <- 1
  result[selective_p_val$change_pts[selective_p_val$pvals <=0.05],6] <- 1

  result[selective_p_val_BH$change_pts[p_selective_p_val_BH <=0.05],7] <- 1
  result[selective_p_val_Bon$change_pts[p_selective_p_val_Bon <=0.05],8] <- 1


  result[,9] <- true_spike
  result[,10] <- as.numeric(dat)
  
  return(result)
}


spike_naive_test <- function(dat, decay_rate, lam, sig, window_size){
  
  fit_spike <- spike_estimates(dat, decay_rate, lam)
  
  spike_loc <- fit_spike$change_pts
  v_collection <- lapply(spike_loc, function(thj)
    construct_v(length(dat), thj, window_size, decay_rate))
  vty <- lapply(v_collection, function(v)sum(v*dat))
  vtv <- lapply(v_collection, function(v)sum(v*v))
  
  p_val_result <- unlist(lapply(seq_along(vty), function(i)naive_z_test(vty[[i]], vtv[[i]],
                                                                        sigma2 = sig)))
  sign_spike <- fit_spike$spike_sign
  
  p_bonferroni <- p.adjust(p_val_result, method="bonferroni")
  p_BH <- p.adjust(p_val_result, method="BH")
  
  return(list(spike_loc = spike_loc,sign_spike=sign_spike,p_val=p_val_result,
              p_bonferroni = p_bonferroni, p_BH = p_BH))
  
}

JNFL_var_estimator <- function(y, decay_rate){
  lag_1_diff <- (y[2:length(y)]-decay_rate*y[1:(length(y)-1)])/sqrt(2)
  lag_2_diff <- (y[3:length(y)]-decay_rate*decay_rate*y[1:(length(y)-2)])/sqrt(2)
  sigma_hat <- sqrt(max(0, 2*var(lag_1_diff)-var(lag_2_diff)))
  if(sigma_hat==0){
    return(MAD_var_estimator(y,decay_rate))
  }
  return(sigma_hat)
}


downsample <- function(y, factor) {
    x <- rep(1, factor)
    out <- convolve(y, rev(x), type = 'filter')
    return(out[seq(1, by = factor, to = length(out))])
}

corr_metric <- function(times_x, times_y, factor) {
    times_x_down <- downsample(times_x, factor)
    times_y_down <- downsample(times_y, factor)
    if (sum(times_x_down) == 0 || sum(times_y_down) == 0) {
        return(0)
    }
    return(cor(times_x_down, times_y_down))
}


cppFunction('double vanRossumDist_cpp(NumericVector u, NumericVector v, double tau) {
            int len_u = u.size();
            int len_v = v.size();
            double total = 0;

            for(int i = 0; i < len_u; ++i) {
            for (int j = 0; j < len_u; ++j){
            total += exp(-abs(u[i]-u[j])/tau);
            }
            }

            for(int i = 0; i < len_u; ++i) {
            for (int j = 0; j < len_v; ++j){
            total -= 2*exp(-abs(u[i]-v[j])/tau);
            }
            }


            for(int i = 0; i < len_v; ++i) {
            for (int j = 0; j < len_v; ++j){
            total += exp(-abs(v[i]-v[j])/tau);
            }
          }

            return total;
        }')




cppFunction('double VictorPurpuraDist_cpp(NumericVector u, NumericVector v, double cost) {

            int len_u = u.size();
            int len_v = v.size();
            double total = 0;
            double choice_1, choice_2, choice_3;

            if (cost == 0){
              if(len_u >= len_v){
                return(len_u-len_v);
              }else{
                return(len_v-len_u);
              }
            }

            NumericMatrix cost_mat(len_u+1, len_v+1);

            for (int i = 0; i<= len_u; i++){
                cost_mat(i,0) = i;
            }
            for (int i = 0; i<= len_v; i++){
                cost_mat(0,i) = i;
            }

            for (int i = 1; i <= len_u; i++){
              for (int j = 1; j <= len_v; j++){
                choice_1 = cost_mat(i-1,j) + 1;
                choice_2 = cost_mat(i,j-1) + 1;
                if (u[i-1]-v[j-1] >= 0){
                  choice_3 = cost_mat(i-1, j-1) + cost*(u[i-1]-v[j-1]);
                }else{
                  choice_3 = cost_mat(i-1, j-1) + cost*(-u[i-1]+v[j-1]);
                }

                cost_mat(i,j) = std::min(choice_1, std::min(choice_2, choice_3));
              }
            }

            total = cost_mat(len_u,len_v);
            return total;
            }')

