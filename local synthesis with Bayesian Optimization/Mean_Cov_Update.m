%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function updates the mean and covariance of a Gaussian process as a function of  
% the scheduling variable of an LPV system. 
%
% Input:   Kstar: covarince between a new grid point and the observation
% points. Mu2: observed values of the unknown function to be regressed,
% Cov_mat: the covarince matrix of the known points.

% output:  Mu_star: updated mean of the process.( parametrized estimation of the function)
%          sigma_star: uncertainty of estimation. 
% Bardia Sharif
% Eindhoven University of Technology
% August 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Mu_star,sigma_star] =    Mean_Cov_Update(Kstar,Mu2,Cov_mat)
Mu_star      =  Kstar'/Cov_mat*(Mu2);
sigma_star   =  -Kstar'/(Cov_mat)*Kstar + 1; 
end