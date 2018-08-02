%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function produces the  the Dual LMI together with the decision varibales for local 
% LPV synthesis with PILF.  
% 
% Inputs : system interconnection (sysic), number of grid points, number of measured outputs,
% number of control inputs
% 
% reference: A Survey of Linear Parameter-Varying Control
% Applications Validated by Experiments or
% High-Fidelity Simulations, Christian Hoffman
%
% Bardia Sharif
% November 2017
% Eindhoven university of technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [LMI_d, Y, gamma] = Dual_LMI_LOCAL(P,np,ny,nu)
%% get the generalized plant matrices
[A,B1,B2,C1,C2,D1,D12,D21,D22,nx,nz,nw] = get_sysic_matrices(P,np,ny,nu);
%% Construct the LMI blocks
%% Construct the middle term of the LMI for each scheduling varibale (based on eq 17, of the reference paper) 
Y                           =       sdpvar(nx,nx);                                                      % storage function in the primal LMI.
gamma                       =       sdpvar( 1, 1);                                                      % induced L2 gain upper bound.
storage_blk_dual            =       [zeros(nx),Y;Y,zeros(nx)];                                          % Storage function block.

SP                          =       [storage_blk_dual,zeros(size(storage_blk_dual,1),nz);...            % Storage+ performance block.
                                    zeros(nz,size(storage_blk_dual,2)),gamma*eye(nz)];
LMI_d                       =       [];                
for i = 1:np                    
    clear blk_11 blk_12 blk_21 blk_22 buffer middle buffer Phi PHI;
    Quad_p(:,:,i)           =       [-A(:,:,i)'                     , -C1(:,:,i)';...
                                    eye(nx,size(-A(:,:,i)',2))      , zeros(nx,size(C1(:,:,i)',2));...
                                    zeros(nz,size(A(:,:,i)',2))     , eye(nz,size(C1(:,:,i)',2))];
    blk_11                  =       Quad_p(:,:,i)'*SP*Quad_p(:,:,i);                                    % The 1,1 block of the middle expression.
    blk_12                  =       -[B1(:,:,i);D1(:,:,i)];                                             % The 1,2 block of the middle expression.   
    blk_21                  =       -[B1(:,:,i)' D1(:,:,i)'];                                           % The 2,1 block of the middle expression.
    blk_22                  =       gamma*eye(size(blk_21,1),size(blk_12,2));                           % The 2,2 block of the middle expression.
    middle                  =       [blk_11,blk_12;blk_21,blk_22];                                      % Construct the middle expression for this specific value of scheduling variable( grid point).
    Phi                     =       null([B2(:,:,i)' D12(:,:,i)']);
    PHI                     =       [Phi,zeros(size(Phi,1),size(blk_21,1));zeros(size(blk_12,2),size(Phi,2)),eye(size(blk_12,2),size(blk_21,1))];        
    elim_d                  =       PHI;
    buffer                  =       elim_d'*middle*elim_d;
    LMI_d                   =       blkdiag(LMI_d,buffer);
end                  
end