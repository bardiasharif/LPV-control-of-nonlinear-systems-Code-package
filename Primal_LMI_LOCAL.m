%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function produces the  the primal LMI together with the decision varibales for local 
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
function [LMI_p, X, gamma] = Primal_LMI_LOCAL(P,np,ny,nu)
%% get the generalized plant matrices
[A,B1,B2,C1,C2,D1,D12,D21,D22,nx,nz,nw] = get_sysic_matrices(P,np,ny,nu);
%% Construct the LMI blocks
%% Construct the middle term of the LMI for each scheduling varibale (based on eq 16, of the reference paper) 
X                           =       sdpvar(nx,nx);                                                      % storage function in the primal LMI.
gamma                       =       sdpvar( 1, 1);                                                      % induced L2 gain upper bound.
storage_blk_primal          =       [zeros(nx),X;X,zeros(nx)];                                          % Storage function block.

SP                          =       [storage_blk_primal,zeros(size(storage_blk_primal,1),nw);...        % Storage+ performance block.
                                    zeros(nw,size(storage_blk_primal,2)),-gamma*eye(nw)];
LMI_p                       =       [];                                
for i = 1:np                    
    clear blk_11 blk_12 blk_21 blk_22 middle Psi PSI;
    Quad_p(:,:,i)           =       [eye(nx,size(A(:,:,i),2))       , zeros(nx,size(B1(:,:,i),2));...   % The Quadratic part of the 1,1 block of the middle experssion
                                    A(:,:,i)                        , B1(:,:,i);...
                                    zeros(nw,size(A(:,:,i),2))      , eye(nw,size(B1(:,:,i),2))];
    blk_11                  =       Quad_p(:,:,i)'*SP*Quad_p(:,:,i);                                    % The 1,1 block of the middle expression.
    blk_12                  =       [C1(:,:,i)';D1(:,:,i)'];                                            % The 1,2 block of the middle expression.   
    blk_21                  =       [C1(:,:,i) D1(:,:,i)];                                              % The 2,1 block of the middle expression.
    blk_22                  =       -gamma*eye(size(blk_21,1),size(blk_12,2));                          % The 2,2 block of the middle expression.
    middle                  =       [blk_11,blk_12;blk_21,blk_22];                                      % Construct the middle expression for this specific value of scheduling variable( grid point).
    Psi                     =       null([C2(:,:,i) D21(:,:,i)]);
    PSI                     =       [Psi,zeros(size(Psi,1),size(blk_21,1));zeros(size(blk_12,2),size(Psi,2)),eye(size(blk_12,2),size(blk_21,1))];        
    elim_p                  =       PSI;
    buffer                  =       elim_p'*middle*elim_p;
    LMI_p                   =       blkdiag(LMI_p,buffer);
end                  
end