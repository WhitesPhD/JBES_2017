## Replication code

This folder contains the baseline code for the paper "Macroeconomic Factors Strike Back: A Bayesian Change Point Model of Time-Varying Risk Exposures and Premia in the U.S. Cross Section" (joint with Massimo Guidolin and Francesco Ravazzolo) published in 2017 in the Journal of Business Economic and Statistics. 

Author: Bianchi Daniele, daniele.bianchi@qmul.ac.uk

------------------------------------------------------------------------------------------------

### Files description:
------------------------------------------------------------------------------------------------

1. Data.mat: Contains the time-series data about both asset/portfolios and macroeconomic risk factors/instruments

2. Main.m is some of the main matlab file reads the data and run the main algorithm

3. StepUnique.m is the matlab function running the Gibbs sampler

4. breaks_sampler_aff.m runs the Gerlach et al. (2000) algorithm to simulate structural breaks in the betas

5. breaks_sampler_aff_R.m runs the Gerlach et al. (2000) algorithm to simulate structureal breaks in the log volatility

6. forward_filtering_backward_sampler.m runs a standard FFBS algorithm to sample from the posterior of betas and log-volatilities

7. mixtures.m auxiliary function to sample from the chi2(1) 

8. randgamma.m auxiliary function to generate N gamma random variables

Notice: to run some of the function you need the Matlab statistics toolbox
