%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function constructs the permuted multiplier blocks required for
% applying the linearization lemma for constructing a global LPV controller
% via the Full block S procedure
%
% reference: LPV control and full block multipliers by Carsten Scherer/
% Linear matrix inequalities in control by Siep Weiland and Carsten
% Scherer.
% 
% Bardia Sharif
% August 2018
% Eindhoven university of technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Qbar,Sbar,Rbar,Rv] = Permutated_Multipliers(Q,S,R,P,Pe,Xx,nzc,nz,nw,gamma)

Q       =   value(Q);
S       =   value(S);
R       =   value(R);
St      =   P(size(Q,1)+1:end,1:size(S,1));

S12     =   Pe(1:nzc,end-nzc+1:end);
Q12     =   Pe(1:size(Q,1),size(Q,2)+size(S,2)+1:end-nzc);

S21t    =   Pe(size(Q,1)+1:size(St,1)+size(Q,1),size(St,2)+size(R,2)+1:size(St,2)+size(R,2)+size(Q12,2));
R12     =   Pe(size(Q,1)+1:size(St,1)+size(Q,1),size(St,2)+size(R,2)+size(S21t,2)+1:size(St,2)+size(R,2)+size(S21t,2)+size(S12,2));

Q12t    =   Pe(size(Q,1)+size(St,1)+1:end-size(S12',1),1:size(Q,2));
S21     =   Pe(size(Q,1)+size(St,1)+1:end-size(S12',1),size(Q,2)+1:size(Q,2)+size(S,2));

Q2      =   Pe(size(Q,1)+size(S',1)+1:end-size(S12',1),size(Q12',2)+size(S21,2)+1:size(Q12',2)+size(S21,2)+size(Q12,2));
S2      =   Pe(size(Q,1)+size(S',1)+1:end-size(S12',1),size(Q12',2)+size(S21,2)+size(Q2,2)+1:end);

S12t    =   Pe(size(Q,1)+size(S',1)+size(Q12',1)+1:end,1:size(Q,2));
R12t    =   Pe(size(Q,1)+size(S',1)+size(Q12',1)+1:end,size(Q,2)+1:size(Q,2)+size(S,2));
S2t     =   Pe(size(Q,1)+size(S',1)+size(Q12',1)+1:end,size(Q,2)+size(R12t,2)+1:size(Q,2)+size(R12t,2)+size(Q12,2));

R2      =   Pe(size(Q,1)+size(S',1)+size(Q12',1)+1:end,size(S12',2)+size(R12',2)+size(S2',2)+1:end);
%%  Construct the Qbar, Sbar, Rbar for applying the linearization lemma
qq      =   [Q Q12;Q12t Q2];
Qbar    =   blkdiag(qq,-gamma*eye(nw));
Qbar    =   blkdiag(zeros(size(Xx,1),size(Xx,2)),Qbar);

ss      =   [S S12;S21 S2];
Sbar    =   blkdiag(ss,zeros(nw,nz));
Sbar    =   blkdiag(Xx,Sbar);

rr      =   [R R12;R12t R2];
Rv      =   blkdiag(rr,(1/gamma)*eye(nz));
Rbar    =   blkdiag(zeros(size(Xx,1),size(Xx,2)),Rv);
end