%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function constructs the LTI block of a global LPV controller by
% applying Linearization lemma to eq22 of the reference paper.
%
% inputs: Matrices of the closed-loop system (eq.18), permuted multipliers 
% obtained from the function "Permutated_Multipliers.m", the first term and left and right
% matrices of the second term of Eq18 from the reference paper which can be
% obtained from the function "Closed_Loop_Matrices.m"
%
% reference: LPV control and Full block multipliers: Carsten Scherer
%
% Bardia Sharif
% August 2018
% Eindhoven university of technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Gk = Controller_Construction(Gcl,Qbar,Sbar,Rbar,Rv,left_mat,right_mat,first_term )
%% construction of matrices required for applying the linearization lemma.
T   =   chol(Rv,'lower');
Uv  =   eye(size(T,2));
%% set up yalmip options
yalmip('clear');
ops          = sdpsettings('savesolveroutput', 1, ...
               'savesolverinput' , 1, ...
               'verbose'         , 1, ...
               'solver'          ,'sdpt3');  
%% Constructed the linearized QMI and solve it to get the closed loop matrices ( Lemma4.1 : Schrer and Weiland)
LMI_eq  =   [Qbar + Sbar*Gcl + Gcl'*Sbar', Gcl'*T;...
             T'*Gcl, -Uv];
cons    =   LMI_eq <= 0;
sol     =   optimize(cons,[],ops);
%% Use the solution of the LMI to construct the LTI part of the controller 
Gk  =   (left_mat)\(Gcl-first_term)/right_mat;
end