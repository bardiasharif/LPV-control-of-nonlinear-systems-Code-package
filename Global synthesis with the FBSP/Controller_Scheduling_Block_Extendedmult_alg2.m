%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function produces the scheduling block Delta_c , the extended multiplier Pe and the 
% closed loop lyapunov finction Xx  of a global LPV
% controller based on Algorithm 2 of the reference paper
% 
% Inputs : Multiplierss, system scheduling block
% 
% reference : LPV control and full block multipliers by Carsten Scherer
%
% Bardia Sharif
% Jannuary 2018
% Eindhoven university of technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
function [Delta_c,U,T,Nmin,Nplus,Vmin,Vplus,Pe,Xx,Kc,mc]    =   Controller_Scheduling_Block_Extendedmult_alg2(P,Ptil,X,Y,DELTA)
%% perturb Ptil if necessary to render it non-singular
Ptil                        =   value(Ptil);

%% Construct the scheduling block of the controller, the extended multiplier and the Lyapunov variable \cal{X} (Naming convention follows from the reerence paper)
P       =   value(P);
X       =   value(X);
Y       =   value(Y);
diff    =   P-inv(Ptil);
T       =   eye(size(P-inv(Ptil)));
Sdelta  =   [DELTA;eye(size(DELTA))];


[U,D]             =   eig(diff);
[~,permutation]   =   sort(diag(D));                    % order the eigenvalues in ascending order.
D                 =   D(permutation,permutation);       % order the eigenvalues in ascending order.   
U                 =   U(:,permutation);                 % reorder the eigenvectors for diagonalizing N in ascending order.
N                 =   U'*diff*U;


tol               =  1e-8;
jj                =  find(abs(N)<tol);
N(jj)             =  0;                                 % Set machine precision zeros to 0.

Kc     =   length(find(N>tol));
mc     =   length(find(N<-tol));

Nmin   =   N(1:mc,1:mc);
Nplus  =   N(mc+1:end,mc+1:end);

V       =   Sdelta'*U;
Vmin    =   V(:,1:mc);

Pe      =   [P U*T;T'*U' T'*(pinv(N))*T];
% [Z,~]   =   eig(X-inv(Y));
Z       =   orth(X-inv(Y));

Vplus  =   V(:,mc+1:end);

Delta_c     =  Nmin*Vmin'*inv(Sdelta'*P*Sdelta - Vmin*Nmin*Vmin')*Vplus;

Xx      =   [X Z;Z' pinv(Z'*(X-inv(Y))*Z)];
end