

packages <- c("devtools", "tidyverse", "latex2exp", "RColorBrewer",
              "shades", "ggpubr")
install.packages(setdiff(packages, rownames(installed.packages())))  
if(!("SpikeInference" %in% rownames(installed.packages()))){
  devtools::install_github("yiqunchen/SpikeInference")
}
