%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function produces the primal LMI for global LPV synthesis based on the reference paper 
% Eq 25. 
% 
% Inputs : Multiplier,system interconnection (sysic), number of control inputs and
% number of measurement outputs
% 
% reference : LPV control and full block multipliers by Carsten Scherer
%
% Bardia Sharif
% Jannuary 2018
% Eindhoven university of technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
function [LMI_Primal,X,gamma]     =   LMI_Primal_Global(Q,S,R,Pw,ny,nu)
%% Get the system matrices based on the convention in the reference
[A,Bu,Bp,B,Cu,Cp,C,Duu,Dup,Eu,Dpu,Dpp,Ep,Fu,Fp,D,nx,nzu,nwu,nw,nz,DELTA] = get_sys_matrices_global(Pw,ny,nu);
%% initialize the LMI blocks
Psi             =   null([C Fu Fp]);    % The basis function for applying the elimination lemma.

Middle_term     =   [];                 % The middle term and its blocks
%% construct the middle term
X               =   sdpvar(nx,nx);                  % The storage function
gamma           =   sdpvar(1,1);                    % Upper bound on the L2 induced gain of the system

blk11   =       [zeros(nx,nx),X; X,zeros(nx,nx)];
blk12   =       [zeros(size(X,1),size(Q,2)),zeros(size(X,1),size(S,2));...
                 zeros(size(X,1),size(Q,2)),zeros(size(X,1),size(S,2))];
blk13   =       [zeros(size(X,1),nw);zeros(size(X,1),nw)];

blk21   =       [zeros(size(Q,1),size(X,2)),zeros(size(Q,1),size(X,2));...
                 zeros(size(R,1),size(X,2)),zeros(size(R,1),size(X,2))];
blk22   =       [Q,S;S',R];
blk23   =       [zeros(size(Q,1),nw);zeros(size(R,1),nw)];

blk31   =       [zeros(nw,nx) zeros(nw,nx)];
blk32   =       [zeros(nw,size(Q,2)),zeros(nw,size(S,2))];
blk33   =       -gamma*eye(nw);


middle_term     =   [blk11 blk12 blk13;blk21 blk22 blk23;blk31 blk32 blk33];


blkq11  =   [eye(nx,size(A,2)),zeros(nx,size(Bu,2)),zeros(nx,size(Bp,2));...
             A, Bu, Bp];
         
blkq21  =   [zeros(size(Q,2),size(Cu,2)),eye(size(Q,2),size(Duu,2)),zeros(size(Q,2),size(Dup,2));...
             Cu, Duu, Dup];
         
blkq31  =   [zeros(nw,size(Cu,2)),zeros(nw,size(Duu,2)),eye(nw,size(Dup,2))];


Quad_term       =   [blkq11;blkq21;blkq31];

BLK11           =   Quad_term'*middle_term*Quad_term;
BLK12           =   [Cp'; Dpu';Dpp'];
BLK21           =   BLK12';
BLK22           =   -gamma*eye(nz);

Middle_term     =   [BLK11 BLK12;BLK21 BLK22];

Elim            =   blkdiag(Psi ,eye(nz));

LMI_Primal      =   Elim'*Middle_term*Elim;
end