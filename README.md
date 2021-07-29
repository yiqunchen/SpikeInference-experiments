# SpikeInference-experiments
Code and instructions to reproduce results and figures from the paper:

Chen YT, Jewell SW, Witten DM. (2021) [Quantifying uncertainty in spikes estimated from calcium imaging data](https://arxiv.org/abs/2103.07818v1
). arXiv:2103.0781 [statME].

### Step-by-step instructions for reproducing figures
#### Note that all the R codes assume that the default working directory is `~/Desktop/SpikeInference-experiments` which might not be the case. Please double check your working directories and modify the code accordingly (e.g., if you are running them on a computing cluster).

1. Download this repository, install relevant packages, and create local directories.
```
git clone https://github.com/yiqunchen/SpikeInference-experiments.git
cd SpikeInference-experiments
Rscript install_packages.R
mkdir -p output_files
mkdir -p plot_output
```
2. Generate datasets for Figures 1 and 3(a)-(b) (Empirical Selective Type I Error).
```
Rscript ./Figure_1_3_Gen_Data.R
```
3. Generate Figures 1 and 3(a)-(b).
```
Rscript ./Figure_1_3.R
```
4. Generate datasets for Figures 3(c)-(d) and 4 (Empirical Power and Confidence Interval).
```
Rscript ./Figure_3_4_Gen_Data.R
```
5. Generate Figures 4 and 5.
```
Rscript ./Figure_3_4.R
```
6. Generate Figure 2.
```
Rscript ./Figure_2.R
```
7. Generate intermediate files used for Figure 6.
```
Rscript ./2d_tuning_param_ell_0.R
Rscript ./cluster_qc_spikefinder.R
Rscript ./refine_spike_2d_param.R
```
8. Generate Figure 6.
```
Rscript ./Figure_6.R
```
9. Generate intermediate files used for Figure 7 (Note: this process might be time-consuming; you can reduce the number of resampling distribution accordingly).
```
Rscript ./perm_test_2d_param.R
```
10. Generate Figure 7.
```
Rscript ./Figure_7.R
```
11. Generate Figure 8.

```
Rscript ./timing_data.R
Rscript ./Figure_8.R
```

12. Finally, Figures 9 and 10 can be reproduced by changing a plotting parameter (*h*) in step 10.


