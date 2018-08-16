%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function produces the clsoed system matrices of an LPV system in the generalized plant
% framework (Based on Eq 18 of the reference paper)as a function of the LTI part of a LPV controller.
% 
% reference : LPV control and full block multipliers by Carsten Scherer
%
% Bardia Sharif
% Jannuary 2018
% Eindhoven university of technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
function [CL_mat,left_mat,right_mat,first_term,nc,nwc,nzc] = Closed_Loop_Matrices(A,Bu,Bp,B,Cu,Cp,C,Duu,Dup,Eu,Dpu,Dpp,Ep,Fu,Fp,D,ny,nu,nx,nzu,nwu,nw,nz,mc,Kc,X,Y)
X   =   value(X);
Y   =   value(Y);

nc      =   length(range(X-inv(Y)));
nwc     =   mc;
nzc     =   Kc;

% parametrize the controller matrices
Ac      =   sdpvar(nc,nc,'full');
Bc1     =   sdpvar(nc,ny,'full');
Bc2     =   sdpvar(nc,nwc,'full');

Cc1     =   sdpvar(nu,nc,'full');
Dc1     =   sdpvar(nu,ny,'full');
Dc12    =   sdpvar(nu,nwc,'full');

Cc2     =   sdpvar(nzc,nc,'full');
Dc21    =   sdpvar(nzc,ny,'full');
Dc2     =   sdpvar(nzc,nwc,'full');

Kmat    =   [Ac,Bc1,Bc2;...
             Cc1,Dc1,Dc12;...
             Cc2,Dc21,Dc2];
%% Write the closed loop matrices as a function of controller matrices (Eq18 Carsten's paper)
left_mat        =   [zeros(size(B,1),nc),B,zeros(size(B,1),Kc);...
                     eye(nc),zeros(nc,size(B,2)),zeros(nc,Kc);...
                     zeros(size(Eu,1),nc),Eu,zeros(size(Eu,1),Kc);...
                     zeros(Kc,nc),zeros(Kc,size(Eu,2)),eye(Kc);...
                     zeros(size(Ep,1),nc),Ep,zeros(size(Ep,1),Kc)];
                 
right_mat       =    [zeros(nc,size(C,2)),eye(nc),zeros(nc,size(Fu,2)),zeros(nc,mc),zeros(nc,size(Fp,2));...
                      C,zeros(size(C,1),nc),Fu,zeros(size(Fu,1),mc),Fp;...
                      zeros(mc,size(C,2)),zeros(mc,nc),zeros(mc,size(Fu,2)),eye(mc),zeros(mc,size(Fp,2))];
  
second_term     =    left_mat*Kmat*right_mat;

ncl             =    nx+nc;

first_term      =    [A, zeros(nx,nc), Bu, zeros(nx,nwc),Bp;...
                      zeros(nc,nx), zeros(nc,nc),zeros(nc,size(Bu,2)),zeros(nc,nwc),zeros(nc,size(Bp,2));...
                      Cu, zeros(size(Cu,1),nc), Duu, zeros(size(Cu,1),nwc), Dup;...
                      zeros(nzc,nx), zeros(nzc,nc), zeros(nzc,size(Duu,2)), zeros(nzc,nwc),zeros(nzc,size(Bp,2));...
                      Cp, zeros(size(Cp,1),nc), Dpu, zeros(size(Cp,1),nwc),Dpp];
 
CL_mat          =    first_term+ second_term;
end