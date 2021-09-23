
------------------------------------------------------------------------------------------------
Paper title: "Macroeconomic factors strike back: A Bayesian change-point model of time-varying 
              risk exposures and premia in the U.S. cross-section"
------------------------------------------------------------------------------------------------

Authors: Bianchi Daniele, daniele.bianchi@unibocconi.it
	 Guidolin Massimo, massimo.guidolin@unibocconi.it
	 Ravazzolo Francesco, Francesco.Ravazzolo@Norges-Bank.no

To report bugs: 

daniele.bianchi@wbs.ac.uk

------------------------------------------------------------------------------------------------

Files description:
------------------------------------------------------------------------------------------------

1. Data.mat: Contains the time-series data about both asset/portfolios and macroeconomic risk factors/instruments

2. Main.m is the main matlab file reads the data and run the main algorithm

3. StepUnique.m is the matlab function running the Gibbs sampler

4. breaks_sampler_aff.m runs the Gerlach et al. (2000) algorithm to simulate structural breaks in the betas

5. breaks_sampler_aff_R.m runs the Gerlach et al. (2000) algorithm to simulate structureal breaks in the log volatility

6. forward_filtering_backward_sampler.m runs a standard FFBS algorithm to sample from the posterior of betas and log-volatilities

7. mixtures.m auxiliary function to sample from the chi2(1) 

8. randgamma.m auxiliary function to generate N gamma random variables

------------------------------------------------------------------------------------------------

Notice: to run some of the function you need the Matlab statistics toolbox
------------------------------------------------------------------------------------------------
