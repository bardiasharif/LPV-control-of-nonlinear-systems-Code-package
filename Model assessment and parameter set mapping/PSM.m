%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used for reducing the parameter set of an affine LPV model 
% based on trajectories experienced during the operation of the system.
% 
% Input: The vecorized state sapce matrices "Gamma" and typical scheduling
% variables "rho_typical"
%
% reference : LPV control of nonlinear systems, Bardia Sharif
%
% Bardia Sharif
% Eindhoven University of Technology
% August 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T,delta]   = PSM(Gamma,rho_typical)

%% perform PCA
[Un,Sn,Vn]  =   svd(Gamma,'econ');
%% Ask the order of trunction from the user based on the singular values
disp(Sn);
prompt = 'Specify how many scheduling variables should be kept? ';
r     =   input(prompt);
%% construct the truncated matrices
Sn    =   Sn(1:r,1:r);
Un    =   Un(:,1:r);
Vn    =   Vn(:,1:r);
%% construct the new scheduling variable and the mapping from the measurable signal and this new parameter
delta       =   Sn*Vn';
T           =   (Sn*Vn'*rho_typical')/(rho_typical*rho_typical');
end