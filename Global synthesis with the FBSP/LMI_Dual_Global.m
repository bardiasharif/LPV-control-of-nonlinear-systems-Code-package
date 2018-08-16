%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function produces the dual LMI for global LPV synthesis based on the reference paper 
% Eq 26. 
% 
% Inputs : Multiplier blocks,system interconnection (sysic), number of control inputs and
% number of measurement outputs
% 
% reference : LPV control and full block multipliers by Carsten Scherer
%
% Bardia Sharif
% Jannuary 2018
% Eindhoven university of technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
function [LMI_Dual,Y]     =   LMI_Dual_Global(Qtil,Stil,Rtil,Pw,ny,nu,gamma)

%% Get the system matrices based on the convention in the reference
[A,Bu,Bp,B,Cu,Cp,C,Duu,Dup,Eu,Dpu,Dpp,Ep,Fu,Fp,D,nx,nzu,nwu,nw,nz,DELTA] = get_sys_matrices_global(Pw,ny,nu);
%% initialize the LMI blocks
Phi             =   null([B' Eu' Ep']);    % The basis function for applying the elimination lemma.

Quad_term       =   [];                    % The quadratic term of the LMI
Middle_term     =   [];                    % The middle term and its blocks
%% construct the middle term of the 11 block
Y               =   sdpvar(nx,nx);                  % The storage function

blk11   =   [zeros(nx,nx),Y;Y,zeros(nx,nx)];
blk12   =   [zeros(nx,size(Qtil,2)),zeros(nx,size(Stil,2));...
             zeros(nx,size(Qtil,2)),zeros(nx,size(Stil,2))];
blk13   =   [zeros(nx,nz);zeros(nx,nz)];


blk21   =   [zeros(size(Qtil,1),nx),zeros(size(Qtil,1),nx);...
             zeros(size(Rtil,1),nx),zeros(size(Rtil,1),nx)];
blk22   =   [Qtil,Stil;Stil',Rtil];
blk23   =   [zeros(size(Qtil,1),nz);zeros(size(Rtil,1),nz)];


blk31   =   [zeros(nz,nx) zeros(nz,nx)];
blk32   =   [zeros(nz,size(Qtil,2)) zeros(nz,size(Stil,2))];
blk33   =   gamma*eye(nz);

middle_term     =   [blk11 blk12 blk13;blk21 blk22 blk23; blk31 blk32 blk33];
%% quadratic term
blkq11  =   [-A' -Cu' -Cp'; eye(nx,size(A',2)) zeros(nx,size(Cu',2)) zeros(nx,size(Cp',2))];
blkq21  =   [-Bu' -Duu' -Dpu';...
             zeros(size(Stil,2),size(Bu',2)), eye(size(Stil,2),size(Duu',2)),zeros(size(Stil,2),size(Dpu',2))];
blkq31  =   [zeros(nz,size(A',2)),zeros(nz,size(Cu',2)),eye(nz,size(Cp',2))];

Quad_term   =   [blkq11;blkq21;blkq31];

%% Construct the middle term of the LMI
BLK11   =   Quad_term'*middle_term*Quad_term;
BLK12   =   [-Bp;-Dup;-Dpp];
BLK21   =   BLK12'; 
BLK22   =   gamma*eye(nw); 
Middle_term     =   [BLK11 BLK12;BLK21 BLK22];
%% Construct the whole LMI
Elim        =   blkdiag(Phi,eye(nw));
LMI_Dual    =   Elim'*Middle_term*Elim;
end