%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function checks quadratic detectability of an affine LPV model at
% the vertices
%
% Input : system interconnection(sysic), number of vertices(nv), number of
% measured outputs(ny), number of control inputs(nu)
% 
% reference paper: Self-scheduled %!& Control of Linear Parameter-varying
% Systems: a Design Example by Pierre Apkarian
%
% Bardia Sharif
% Eindhoven University of Technology
% March 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function detectcheck = quadratic_detectability_check(P,nv,ny,nu)

[A,B1,B2,C1,C2,D1,D12,D21,D22,nx,nz,nw] = get_sysic_matrices(P,nv,ny,nu);

%% declare the necessary variables and matrices
yalmip('clear');
ops          = sdpsettings('savesolveroutput', 1, ...
               'savesolverinput' , 1, ...
               'verbose'         , 1, ...
               'solver'          ,'sdpt3');
X       =  sdpvar(nx);              %declare the lyapunov varaible.
cons    =  X>=0;  
%% construct the LMI constraints
for i = 1:nv
   N            =   null(C2(:,:,i));
   LMI_stab     =   N'*(A(:,:,i)*X + X*A(:,:,i)')*N;
   cons         =  [cons,LMI_stab<=0];       
end
%% solve the LMI
sol                 =   optimize(cons,[],ops);
%% if the LMI is succesfully solved return 1 indicating detectability otherwise return 0
if (sol.problem==0)
    detectcheck   =   1;
else
    detectcheck =     0;
end

end