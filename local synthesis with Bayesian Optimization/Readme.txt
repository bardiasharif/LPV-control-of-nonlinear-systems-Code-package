This folder contains four functions used for grid point selection by means of Bayesian optimization.


1- "log_likelihood.m":
	
	This function constructs the logarithm likelihood function as a function of hyper-parameters. By maximizing 
	this function over the hyper-parameters, a good estimation of the hyper-parameters can be obtained.

2- "squared_exponential_kernel.m":
	
	This function uses the squared exponential kernel for construction of the covariance matrix of a Gaussian
	process. 

3- "Augment_Covariance.m":
 	
	This function is used for augmenting the covariance matrix ( Obtained from squared_exponential_kernel.m) as a 
	function of new grid point.

4- "Mean_Cov_Update.m":
	
	This function uses the result of "Augment_Covariance.m" to update the mean and covariance of the Gaussian process
	as a function of the new parameter grid point. 

Once an estimation of the mean and variance of the process is obtained, the next grid point can be chosen as one which
maximizes a desired utility function.