-This folders contains four functions:

1-  "get_sysic_matrices.m":
 
This function is used for finding the state space matrices of a system interconnection (sysic object in MATLAB)
in the generalized plant framework for a gridded LPV model. It is necessary to find these matrices in order to construct the 
LMIs for controller synthesis.


2-   "Primal_LMI_LOCAL.m":

This function is used for construction of the primal LMI for local (gridded) LPV synthesis.

** reference: A Survey of Linear Parameter-Varying Control Applications Validated by Experiments or High-Fidelity Simulations, 
   Christian Hoffman.


3-   "Dual_LMI_LOCAL.m":

This function is used for construction of the dual LMI for local (gridded) LPV synthesis.

** reference: A Survey of Linear Parameter-Varying Control Applications Validated by Experiments or High-Fidelity Simulations, 
   Christian Hoffman.

4-    "Construct_Controller_local.m":

This function uses the outputs of the rest of the functions in this folder in order to contruct a gridded LPV controller

** reference: Advanced Gain-Scheduling Techniques for Uncertain Systems,Pierre Apkarian.