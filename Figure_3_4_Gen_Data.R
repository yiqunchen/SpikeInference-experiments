# some simple power analysis
# parameters to vary
# h: width of the window
# mu: poisson means
# need to tune lambda so that we have abt the right # of cps

###### no cp at all ######
library(SpikeInference)

helper_function_dir <- "~/Desktop/SpikeInference-experiments/"
output_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
input_dir <- "~/Desktop/SpikeInference-experiments/input_files/"
plot_output_dir <- "~/Desktop/SpikeInference-experiments/plot_output/"

source(paste0(helper_function_dir,"spike_inf_helper.R"))

set.seed(1995)
# parameters to vary
#change sim_times to 500 to reproduce paper figures
sim_times <- 10
gam <- 0.98
h_tol <- 2
alpha_level <- 0.05
n_length <- 10000
sigma_sim <- rev(round(1/c(1:10),3))

firing_rate_vec <- c(0.01)
h_vec <- c(1,2,10,20)

for (firing_rate in firing_rate_vec){
  target_rate <- firing_rate

  for (h in h_vec){
    for (sigma in sigma_sim){

      sign_vec <- vector("list", length = sim_times)
      pval_vec <- vector("list", length = sim_times)
      j_i_vec <- vector("list", length = sim_times)
      rej_vec <- vector("list", length = sim_times)
      lcb_vec <- vector("list", length = sim_times)
      ucb_vec <- vector("list", length = sim_times)
      true_spike_vec <- vector("list", length = sim_times)
      exp_id <- vector("list", length = sim_times)
      detect_vec <- vector("list", length = sim_times)
      vTc_vec <- vector("list", length = sim_times)
      vTy_vec <- vector("list", length = sim_times)
      naive_lcb <- vector("list", length = sim_times)
      naive_ucb <- vector("list", length = sim_times)

      for (i in c(1:sim_times)){

        cat('currently at sim ', i, " mean ",firing_rate," sigma ", sigma, "h", h,'\n')

        curr_sim <- simulate_ar1(n = n_length, gam = gam, poisMean = firing_rate,
                                 sd = sigma, seed = i)

        if(sigma>=0.5){
          lam_max = 5
        }else{
          lam_max = 1
        }

        fit_spike <- estimate_spike_by_firing_rate(curr_sim,
                           decay_rate = gam, target_firing_rate = target_rate,
                           lam_min = 1e-7, lam_max = lam_max, max_iters=10,
                           tolerance=max(5,floor(n_length*firing_rate*0.05)))


        inference_spike <- spike_inference(dat = curr_sim$fl, decay_rate = gam,
                                           tuning_parameter = fit_spike$tuning_parameter,
                                           window_size = h,
                                           sig2 = sigma*sigma,
                                           return_conditioning_sets = FALSE,
                                           return_ci = TRUE)

        spike_pvals <- inference_spike$pvals
        names(spike_pvals) <- inference_spike$change_pts

        j_i <- nearest_changepoint(curr_sim$spikes,inference_spike$change_pts)
        rej_criterion_1 <- abs(j_i-curr_sim$spikes)<=h_tol
        rej_criterion_2 <- spike_pvals[as.character(j_i)]<=alpha_level
        rej_test <- (rej_criterion_1&rej_criterion_2)


        pval_vec[[i]] <- spike_pvals
        rej_vec[[i]] <- rep(mean(rej_test),times=length(inference_spike$change_pts))
        detect_vec[[i]] <- rep(mean(rej_criterion_1),times=length(inference_spike$change_pts))
        lcb_vec[[i]] <- inference_spike$LCB
        ucb_vec[[i]] <- inference_spike$UCB
        exp_id[[i]] <- rep(i, times=length(inference_spike$change_pts))

        # compute the collection of v
        v_collection <- lapply(inference_spike$change_pts, function(thj)
          construct_v(length(curr_sim$fl), thj, h, gam))
        vTc_collection <- lapply(v_collection, function(v) sum(v*curr_sim$conc)) # the true parameter is vTc
        true_spike_vec[[i]] <- curr_sim$conc[(inference_spike$change_pts)+1]-gam*curr_sim$conc[(inference_spike$change_pts)]
        vTc_vec[[i]] <- unlist(vTc_collection)

        vTy_vec[[i]] <- unlist(lapply(v_collection, function(v) sum(v*curr_sim$fl)))
        v_norm <- unlist(lapply(v_collection,function(v)sum(v*v)))
        naive_lcb[[i]] <- vTy_vec[[i]]-qnorm(0.975)*sqrt(v_norm*sigma*sigma)
        naive_ucb[[i]] <- vTy_vec[[i]]+qnorm(0.975)*sqrt(v_norm*sigma*sigma)



      }

      df_analysis <- data.frame(
                                  pval_vec=unlist(pval_vec),
                                  sigma_sim = sigma,
                                  detection_prob = unlist(detect_vec),
                                  rejection_pow = (unlist(rej_vec)),
                                  lcb_vec = unlist(lcb_vec),
                                  ucb_vec = unlist(ucb_vec),
                                  vTc = unlist(vTc_vec),
                                  true_spike_vec = unlist(true_spike_vec),
                                  exp_id = unlist(exp_id),
  vTy_vec = unlist(vTy_vec),
                                  naive_lcb = unlist(naive_lcb),
  naive_ucb = unlist(naive_ucb))


      save(df_analysis,
           file = paste0(output_dir,
    "power","_mean_",firing_rate,"_sigma_",sigma,"_h_",h,".RData"))

    }
  }
}














