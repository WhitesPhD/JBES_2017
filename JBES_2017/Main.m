%--------------------------------------------------------------------------
% Main code to run the Bayesian change-point model with time-varying
% betas and conditional volatility
% Refer to it as:
% Bianchi, D., M. Guidolin, and F. Ravazzolo (2013), "Macroeconomic Factors Strike Back: 
% A Bayesian Change-Point Model of Time-Varying Risk Exposures and Premia
% in the U.S. Cross-Section", Norges Bank working paper.
%--------------------------------------------------------------------------
% To report bugs:
% daniele.bianchi@unibocconi.it
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Housekeeping
%--------------------------------------------------------------------------

close all; clear all; clc; pause(0.01), randn('seed',sum(clock)), rand('seed',sum(clock))

%--------------------------------------------------------------------------
% Adding folders 
%--------------------------------------------------------------------------

addpath([pwd '\Data\']);
addpath([pwd '\Functions\']);


%--------------------------------------------------------------------------
% Load dataset: The data on both risk factors and asset/portfolios are
% separated between training (prior calibration) and testing sample.
% 
% DataPrior: Data to calibrated the priors
% DataUse  : Data to feed the model
%--------------------------------------------------------------------------

load Data.mat

%--------------------------------------------------------------------------
% Run the main function: Gibbs Sampler
%--------------------------------------------------------------------------

[MCMC_B, MCMC_BKF, MCMC_K, MCMC_L, MCMC_Q, MCMC_R, MCMC_prob] = StepUnique(DataPrior,DataUse,Hyp,MCMC);

%--------------------------------------------------------------------------
% Outputs:
%--------------------------------------------------------------------------
% MCMC_B   : (M x T x N x K) Array of time-varying betas simulated from the posterior betas
% MCMC_BKF : (M x T x N x K) Array of time-varying betas from the filtering distribution
% MCMC_K   : (M x T x N x K1) Array of posterior breaks for both betas and log volatility
% MCMC_Q   : (M x N x K1) Array of samples from the posterior of breaks size
% MCMC_R   : (T x N x M) Array of samples of the conditional idiosyncratic volatilities
% MCMC_prob: (M x N x K) Array of posterior breaks probabilities for the betas



