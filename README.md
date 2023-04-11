<!--- left_trunc_DR --->

# Doubly Robust Estimation under Covariate-Induced Dependent Left Truncation

This repository contains the code that implements the estimators proposed in the paper *Doubly Robust Estimation under Covariate-Induced Dependent Left Truncation* by Yuyao Wang, Andrew Ying, and Ronghui Xu. 

The folder 'src' contains the functions that implement the estimators. The R scripts under the root directory of this repository contain the main code for simulations and data analysis.
In addition, the folder 'Scenarios1-7' contains code for additional simulations under the 7 scenarios in Section 6.1 of the arXiv preprint arXiv:2208.06836.




## Simulation
 
### Simulation without right censoring
 
 - 'simu_save.R': The code that simulates and saves the data sets.

 - 'main.cox_naive_full.R': The main code that load the simulated data sets and computes the 'dr', 'IPW.Q', 'Reg.T1', 'Reg.T2' estimators with Cox models to estimate the two nuisance parameters, the naive estimator that ignores left truncation, and the oracle full data estimator.
 
 - 'main.ltrcrrf.R': The main code that load the simulated data sets and computes the 'cf' estimator where the nuisance parameters $F$ and $G$ are estimated by the relative risk forest (Yao et. al. 2020). In particular, $F$ and $G$ are estimated using the `ltrcrrf()` function of the 'LTRCforests' R package.
 
 - 'main_pl.R': The main code that load the simulated data sets and computes the product-limit (PL) estimator.
 
 - 'main_visulize_simu_result.R': The main code that visualizes the simulation results.
 


 
### Simulation under right censoring

 - 'generate_seeds.R': Generating random seeds for simulating data sets and saving the generated seeds.

 - 'c1_X.main_simu': Simulation for the 'censoring before truncation' scenario ('c1') where censoring can happen before left truncation. The 'dr', 'IPW.Q', 'Reg.T1', and 'Reg.T2' estimators with Cox models to estimate the nuisance parameters, the 'cf-RF-RF' estimator, and the naive estimator and the oracle full data estimator are computed.
 
 - 'c2.main_simu': Simulation for the 'censoring after truncation' scenario ('c2') where censoring is always after left truncation. The 'dr', 'IPW.Q', 'Reg.T1', and 'Reg.T2' estimators with Cox models to estimate the nuisance parameters, the naive estimator, and the oracle full data estimator are computed.
 
 - The folder 'c1.OSG' contains the code for computing the bootstrap SE for the 'cf-RF-RF' estimator using parallel computing on OSG. 
 



### Code for additional simulations with semiparametric models 

As we mentioned before, the folder 'Scenarios1-7' contains the additional simulations in Section 6.1 of the arXiv preprint arXiv:2208.06836.

 - 'simu_save.R': The code that simulates and saves the data sets.
 
 - 'main.cox_naive_full': The code that loads the simulated data sets and computes the 'dr', 'IPW.Q', 'Reg.T1', and 'Reg.T2' estimators with Cox models to estimate the nuisance parameters, the naive estimator, and the oracle full data estimator.
 
 - 'plot_cox_naive': The code for visualizing the simulation results. 
 

 
 
## Applications
 
 - 'c1.main_CNS.R': The main code for analyzing the CNS lymphoma data.
 
 - Folder 'CNS_data' contains the CNS lymphoma data ("CNS_data.txt") from a study of MTX-based chemotherapy (Wang et al., 2015) and the code for the analysis in Vakulenko-Lagun et. al. (2021). The data and the analysis are available in the Supplementary material of Vakulenko-Lagun et. al. (2021).
 
 - 'c2.main_HAAS': The main code for analyzing the HAAS data. 








