This is all the source code needed to reproduce the results repoted in the paper "Overcoming separation and singular random effects covariance matrices in binomial mixed-effects models", by Tina Košuta, Georg Heinze, Ioannis Kosmidis, Philipp Sterzinger, and Rok Blagus. 

All the code was written by Tina Košuta and Rok Blagus. In case of questions or comments please contact rok.blagus@mf.uni-lj.si. There are no restricutions for using and/or adapting this code and or using/sharing the real data used to illustrate the methods. NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE!  Provided "as is".

This folder contains the following data and files that can be used to reproduce all analyses and figures/tables of the manuscript/Supplementary Material.

It contains two subfolders:

./real_data/: [folder for reproducing the results reported in Section "Motivating example", "Illustration"  and Supplementary Material]
script_real_package.R: running this script performs the analysis for the real data example. The script also produces the figures reported in the manuscript/Supplementary Material. Running the script requires the following R packages: AUGglmmTMB,MASS, brglm2, optimx, numDeriv, glmmTMB, lme4, xtable, blme, mvtnorm, dplyr, ggplot2, tidyr, reformulas.  
mv_MSPAL_illustration_2.R: the source code to fit SPM method, obtained from gitHub (https://github.com/psterzinger/softpen_supplementary), adapted by Rok Blagus to allow the calculation of SEs.
AUGglmmTMB_0.1.0.tar.gz: the R package AUGglmmTMB that implements MPL. The package also includes the real data.
figure_real.pdf: coefficient paths (supplementary material)
figure_real1.pdf: coefficient paths (main document)
lims.pdf: conditional log-likelihood by tau (supplementary material)
lims1.pdf: conditional log-likelihood by tau (main document)

./simulated_data/: [folder for reproducing the results reported in Section "Simulation study" and Supplementary Material]
make_scripts.R: running this scipt produces script_i.R, i in 1,...,80 and scriptses_i.R
read_res.R: running this script reads all results (the txt files resulting from runing script_i.R and scriptses_i.R, i in 1,...,80), produces the figures and tables reported in the manuscript and Supplementary Material. Running this script requires the following R packages: ggplot2 and gridExtra.
source_sim.R: this script is used by script_i.R to perform 1000 simulation runs for the specified simulation scenario: this is the core function to reproduce the simulation study. Here, the data according to the specified scenario are created and different methods are used to estimate the model's parameters. The results for each simulation scenarion and method are saved as txt files, where rows correspond to different simulation runs, and columns to the estimated parameters. Running this script requires the following R packages: MASS, brglm2, optimx, lme4, numDeriv, glmmTMB, mvtnorm, and blme. We performed this part of the analysis in the cluster of CentOS based containers
source_sim_getSe.R: this script is used by scriptses_i.R to perform 1000 simulation runs for the specified simulation scenario. These results are required to evaluate if the estimated fixed effects are on the boundary, as suggested by Geroldinger et al. (see the manuscript).

./simulated_data/figs: [folder that contains all figures, in pdf format, shown the manuscript and Supplementary Material]

./simulated_data/results: [folder that contains all txt files resulting from runing script_i.R and scriptses_i.R, i in 1,...,80]

session info for real data illustration (performed on PC):
R version 4.5.2 (2025-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=English_Slovenia.utf8  LC_CTYPE=English_Slovenia.utf8    LC_MONETARY=English_Slovenia.utf8 LC_NUMERIC=C                     
[5] LC_TIME=English_Slovenia.utf8    

time zone: Europe/Ljubljana
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] dplyr_1.1.4         tidyr_1.3.1         ggplot2_4.0.1       car_3.1-3           carData_3.0-5       xtable_1.8-4        mvtnorm_1.3-3      
 [8] blme_1.0-6          numDeriv_2016.8-1.1 lme4_1.1-38         Matrix_1.7-4        optimx_2025-4.9     MASS_7.3-65         glmmTMB_1.1.14     
[15] AUGglmmTMB_0.1.0   

loaded via a namespace (and not attached):
 [1] generics_0.1.4     sandwich_3.1-1     enrichwith_0.4.0   lattice_0.22-7     magrittr_2.0.4     pracma_2.4.6       grid_4.5.2        
 [8] RColorBrewer_1.1-3 nleqslv_3.3.5      brglm2_1.0.1       nnet_7.3-20        Formula_1.2-5      mgcv_1.9-3         purrr_1.2.0       
[15] scales_1.4.0       abind_1.4-8        reformulas_0.4.4   Rdpack_2.6.4       cli_3.6.5          rlang_1.1.6        rbibutils_2.4     
[22] splines_4.5.2      withr_3.0.2        tools_4.5.2        nloptr_2.2.1       minqa_1.2.8        boot_1.3-32        vctrs_0.6.5       
[29] R6_2.6.1           zoo_1.8-14         lifecycle_1.0.4    pkgconfig_2.0.3    pillar_1.11.1      gtable_0.3.6       glue_1.8.0        
[36] Rcpp_1.1.0         statmod_1.5.1      tidyselect_1.2.1   tibble_3.3.0       farver_2.1.2       nlme_3.1-168       labeling_0.4.3    
[43] TMB_1.9.18         compiler_4.5.2     S7_0.2.1        

session info for simulation study (performed on the cluster of CentOS based containers):
R version 3.6.0 (2019-04-26)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

locale:
[1] C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] blme_1.0-6     lme4_1.1-34    Matrix_1.6-2   mvtnorm_1.0-11 glmmTMB_1.1.11

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.11         lattice_0.20-38     reformulas_0.4.0
 [4] rbibutils_2.3       MASS_7.3-51.4       grid_3.6.0
 [7] nlme_3.1-163        Rdpack_2.6.4        minqa_1.2.5
[10] nloptr_1.2.2.1      boot_1.3-28.1       splines_3.6.0
[13] tools_3.6.0         TMB_1.9.6           numDeriv_2016.8-1.1
[16] compiler_3.6.0      mgcv_1.9-0

