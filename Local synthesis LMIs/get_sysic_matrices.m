%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function produces the system matrices of an LPV system in the generalized plant
% framework given the system interconnection(sysic), the number of grid points,
% number of measured outputs and number of control inputs.
% 
% Inputs : system interconnection (sysic), number of grid points, number of measured outputs,
% number of control inputs
% 
% Bardia Sharif
% November 2017
% Eindhoven university of technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
function [A,B1,B2,C1,C2,D1,D12,D21,D22,nx,nz,nw] = get_sysic_matrices(P,np,ny,nu)
%%
[a,b,c,d]   =   ssdata(P);
nx          =   size(P.a,1);               % # of states
nz          =   size(P.c(:,:,1),1) - ny;   % # of generalized performance signals
nw          =   size(P.b(:,:,1),2) - nu;   % # of generalized disturbance signals
%%
% the generalized plant matrices (based on eq 4.1 phase 1 report)
for     i = 1:np  
    A(:,:,i)     =  a.Data(:,:,i);
    B1(:,:,i)    =  b.Data(1:nx,1:nw,i);
    B2(:,:,i)    =  b.Data(1:nx,nw+1:end,i);
    C1(:,:,i)    =  c.Data(1:nz,1:nx,i);
    C2(:,:,i)    =  c.Data(nz+1:end,1:nx,i);
    D1(:,:,i)    =  d.Data(1:nz,1:nw,i);
    D12(:,:,i)   =  d.Data(1:nz,nw+1:end,i);
    D21(:,:,i)   =  d.Data(nz+1:end,1:nw,i);
    D22(:,:,i)   =  d.Data(nz+1:end,nw+1:end,i);
end
end