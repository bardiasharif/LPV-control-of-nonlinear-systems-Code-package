%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function constructs the augmented covariance matrix for a multivariate Gaussian process
% as a function of grid points using the squared exponential kernel
%
% Input:   grid points
% output:  the covariance matrix
%
% Bardia Sharif
% Eindhoven University of Technology
% February 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,Kstar] = Augment_Covariance(Cov_mat,Gpoints,l_f,sig_f)

% syms grid_star
grid_star   =   sym('grid_star',[size(Gpoints,1),1]);
% grid_star = sdpvar(1,1);
n   =   size(Cov_mat,2);
m   =   n;

for i = 1:m
Kstar(i,:)    =  exp((-1/(2*l_f^2))*(norm(Gpoints(:,i)-grid_star,2)^2));

end
K = [Cov_mat+sig_f*eye(n),Kstar;Kstar',1];
end

