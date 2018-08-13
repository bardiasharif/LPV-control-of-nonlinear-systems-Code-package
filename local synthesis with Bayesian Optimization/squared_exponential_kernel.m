%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function constructs the covariance matrix for a Gaussian process
% using the squared exponential kernel
%
% Input:   grid points
% output:  the covariance matrix and the kernel parameters.
%
% Bardia Sharif
% Eindhoven University of Technology
% February 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = squared_exponential_kernel(Grid_points,l_f)

n = length(Grid_points);    % get the number of grid points.

% Construct the covariance matrix as a function of kernel parameters
for    i = 1:n
for    j = 1:n

    K(i,j)  =  exp((-1/(2*l_f^2))*(norm(Grid_points(:,i)-Grid_points(:,j),2))^2);     % compute the elements of the covariance matrix using the squared exponential kernel.
end
end


end 