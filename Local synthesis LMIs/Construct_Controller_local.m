%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function outputs a local LPV controller in state space form by using
% closed form controller equations.
%
% Inputs : Weighted system interconnection (sysic), number of grid points np, number of measured outputs ny,
% number of control inputs nu ,upper bound on the induced L2 gain of the system gamma and the storage functions . 
% X and Y.
%
% referebce: Advanced Gain-Scheduling Techniques
% for Uncertain Systems by Pierre Apkarian
% 
% Bardia Sharif
% November 2017
% Eindhoven university of technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
function [Kp] = Construct_Controller_local(Pw,np,ny,nu,gamma,X,Y)
%% determine the number of states, generalized performance signals and generalized disturbance signals
nx          =   size(Pw.a,1);               % # of states
nz          =   size(Pw.c(:,:,1),1) - ny;   % # of generalized performance signals
nw          =   size(Pw.b(:,:,1),2) - nu;   % # of generalized disturbance signals
%% get the parameter dependednt system interconnection matrices in eq 2 of the reference paper
 A          =   Pw.A;
 B1         =   Pw.B(1:nx,1:nw);
 B2         =   Pw.B(1:nx,nw+1:end);
 C1         =   Pw.C(1:nz,1:nx);
 D11        =   Pw.D(1:nz,1:nw);
 D12        =   Pw.D(1:nz,nw+1:end);
 C2         =   Pw.C(nz+1:end,1:nx);
 D21        =   Pw.D(nz+1:end,1:nw);
 D22        =   Pw.D(nz+1:end,nw+1:end);
%% Construct the matrices   in eq 15-18 of the reference paper
Dk          =   zeros(size(D12,2),size(D21,1));
Dcl         =   D11 + D12*Dk*D21;

Bmat        =   [zeros(size(D21,1),size(D21',2)),D21,zeros(size(D21,1),size(Dcl',2));...
                D21', -gamma*eye(size(D21',1),size(Dcl,2)),Dcl';...
                zeros(size(Dcl,1),size(D21',2)),Dcl,-gamma*eye(size(Dcl',1),size(Dcl,2))]\(-[C2;B1'*X;C1+D12*Dk*C2]);
Bk_hat      =   Bmat(1:ny,1:nx)';

Cmat        =   [zeros(size(D12',1),size(D12,2)),D12',zeros(size(D12',1),size(Dcl,2));...
                D12,-gamma*eye(size(D12,1),size(Dcl,2)),Dcl;...
                zeros(size(Dcl',1),size(D12,2)),Dcl',-gamma*eye(size(Dcl',1),size(Dcl,2))]\(-[B2';C1*Y;(B1+B2*Dk*D21)']);
Ck_hat      =   Cmat(1:nu,1:nx);

Ak_hat      =   -(A+B2*Dk*C2)' + [X*B1+Bk_hat*D21 (C1+D12*Dk*C2)']*([-gamma*eye(size(Dcl',1),size(Dcl,2)), Dcl';Dcl,-gamma*eye(size(Dcl,1),size(Dcl',2))]\[(B1+B2*Dk*D21)';C1*Y+D12*Ck_hat]);
%% Determine M and N for the factorization problem I-XY = NM'
Nc          =    eye(nx);
Mc          =   (eye(nx) - X*Y)';
%% Construct the controller matrices using eq 9-11 in the reference paper
Ak          =   Nc\(Ak_hat - X*(A-B2*Dk*C2)*Y - Bk_hat*C2*Y - X*B2*Ck_hat)/Mc';
Bk          =   Nc\(Bk_hat - X*B2*Dk);
Ck          =   (Ck_hat - Dk*C2*Y)/Mc';
Kp          =   ss(Ak,Bk,Ck,Dk);
end