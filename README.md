<!--- left_trunc_DR --->

# Doubly Robust Estimation under Covariate-induced Dependent Left Truncation

This repository contains code that implement the estimators proposed in the paper 'Doubly Robust Estimation under Covariate-induced Dependent Left Truncation' by Yuyao Wang, Andrew Ying, and Ronghui Xu. 

The 'src' folder contains the functions that implement the estimators, and the R scripts under the root directory of this repository contains the main code for the simlation and data analysis involved in the paper with details listed below. 

In addition, the 'arXiv_20220814' folder includes the code for the simulation involved in the arXiv version of this paper: arXiv:2208.06836, the details of which are also described below.




## Simulation
 
### Simulatiom without right censoring

 - 'main.cox_naive_full.R': The main code for loading the simulated data sets and compute the dr, IPW.Q, Reg.T1, Reg.T2 estimators with `coxph()` for estimating the two nuisance parameters, as well as the naive estimator and the full estimator.
 
 - 'main.ltrcrrf.R': The main code for loading data and compute the 'cf' estimator with `LTRCforests::ltrcrrf()` for estimating the nuisance parameters
 
 - 'main_pl.R': The main function for computing the product-limit estimator.
 
 - 'main_visulize_simu_result.R': The main code for visualizing the simulation results.
 


 
### Simulation under right censoring

 - 'generate_seeds.R': Generate random seeds for simulating data sets and save the seeds.

 - 'c1.main_simu': The main code for simulation under the censoring scenario ('c1') that censoring can happen before left truncation. The estimators computed are the dr, IPW.Q, Reg.T1, and Reg.T2 estimators with coxph() for estimating the nuisance parameters, the cf-RF-RF estimator, as well as the naive and full estimators. 
 
 - 'c2.main_simu': The main code for simulation under the censoring scenario ('c2') that censoring is always after left truncation. The estimators computed are the dr, IPW.Q, Reg.T1, and Reg.T2 estimators with coxph() for estimating the nuisance parameters, as well as the naive and full estimators. 
 
 
 
 
 
## Applications
 
 - 'c1.main_CNS.R': The main code for analyzing the CNS lymphoma data.
 
 - Folder 'CNS_data': the CNS lymphoma data ("CNS_data.txt") from a study of methotrexate-based chemotherapy (Wang et al., 2015) and the code for analyzing this data set from Vakulenko-Lagun et. al. (2021), obtained from the Supplementary material of Vakulenko-Lagun et. al. (2021).
 
 - 'c2.main_HAAS': The main code for analyzing the HAAS data. 






## Code for the arXiv version of this paper: arXiv:2208.06836

As we mentioned before, the folder `arXiv_20220814' contains the simulation code for the arXiv paper: arXiv:2208.06836.

 - 'simu_save.R' contains the code that simulates and saves the datasets.
 
 - 'main.cox_naive_full' contains the code that loads the simulated data sets, computes the estimates from the proposed estimators and their competitors, and saves the results.
 
 - 'plot_cox_naive' contains the code for visualizing the simulation results. 
 

