%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function produces the multipliers and their corresponding matrix
% inequalities required for global LPV synthesis via the full block S
% procedure(based on Eq23 of the reference paper).
%
% Inputs : Generators of the parameter set polytope, the uncertainty block. 
% 
% 
% reference : LPV control and full block multipliers by Carsten Scherer
%
% Bardia Sharif
% Jannuary 2018
% Eindhoven university of technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
function [Multiplier_constraints,P,Ptil,Q,S,R,Qtil,Stil,Rtil] = Multiplier_constraints_FBSP(Vert,DELTA,epsilon)
%% Construct the multipliers  
Deltavert   =   Vert;
del         =   DELTA.Data.NominalValue;
nd          =   size(del,1);
md          =   size(del,2);                                                   

Q       =   sdpvar(nd,nd,'symmetric');
S       =   sdpvar(nd,md,'full');
R       =   sdpvar(md,md,'symmetric');

P       =   [Q,S;S',R];

Qtil    =   sdpvar(nd,nd,'symmetric');
Stil    =   sdpvar(nd,nd,'full');
Rtil    =   sdpvar(nd,nd,'symmetric');

Ptil    =   [Qtil,Stil;Stil',Rtil];

Multiplier_constraints = [Q<=-epsilon*eye(size(Q)),Rtil>=epsilon*eye(size(Rtil))];

%% Construct the constraints at each vertex and concatenate 

numVert       =   size(Deltavert,2);                                        % number of vertices.

Dvertex       =   zeros(nd);                                                % vertex mask.


for i = 1:numVert
   
    for  j = 1:nd
    Dvertex(j,j)    =   Deltavert(j,i);                                     % Populate the vertex mask diagonally with the generators
    end
    
cons1   =   [Dvertex;eye(md)]'*P*[Dvertex;eye(md)];                         % constraints for P
cons2   =   [eye(nd);-Dvertex']'*Ptil*[eye(nd);-Dvertex'];                  % constraints for Ptilde


Multiplier_constraints     =  [Multiplier_constraints,cons1>=eye(size(cons1))*epsilon,cons2<=eye(size(cons2))*-epsilon];    % Concatenate the constraints.
end

end
