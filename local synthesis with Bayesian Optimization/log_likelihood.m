%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.
% This function constructs the loglikelihood of a gaussion kernel as a
% function of its hyper parameters
%
% Bardia Sharif
% Eindhoven University of Technology
% March 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  L = log_likelihood(X,Grid_points,Mu)
n = length(Grid_points);    % get the number of grid points.

for    i = 1:n
for    j = 1:n

    K(i,j)  =  exp((-1/(2*X(1)^2))*(norm(Grid_points(:,i)-Grid_points(:,j),2))^2);     % compute the elements of the covariance matrix using the squared exponential kernel.
end
end

L =   -0.5*Mu'/(K+X(2)*eye(n))*Mu - 0.5*log(det(K+X(2)*eye(n))) -(n/2)*log(2*pi);


end