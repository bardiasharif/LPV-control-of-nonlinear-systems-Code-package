-This folders contains the following  functions all based on:


LPV control and full block multipliers: Caresten Scherer




1-  "LMI_Primal_Global.m":
 
	This function is used for construction of schured version of primal global analysis LMIs as presented in equation 25 of the reference paper. 


2-   "LMI_Dual_GlobalL.m":

	This function is used for construction of a schured version of dual global analysis LMIs as presented in equation 26 of the reference paper. 


3- "Multiplier_constraints_FBSP.m":
	
	This function is used for construction of multiplier constraints.


4- "Controller_Scheduling_Block_Extendedmult_alg2.m":
	
	This function is used for construction of the scheduling block of the controller as well as the closed loop Lyapunov function and the extended multiplier 
	as presented in Algorithm2 of the reference paper.


5- "Closed_Loop_Matrices.m":

	This function constructs the closed loop matrices as a function of the controller matrices as presented in equation 18 of the reference paper.


6- "Closed_Loop_Matrices.m":
	
	This function is used for permutation of multiplier blocks such that linearization lemma can be applied to the QMI in equation 22.


7- "Controller_Construction.m":
	
	This function applies linearization lemma to the QMI in equation 22 and constructs the LTI block of the controller.