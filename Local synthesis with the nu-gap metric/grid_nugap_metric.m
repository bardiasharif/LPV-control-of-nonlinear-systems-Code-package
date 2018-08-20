%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function grids the parameter space of an LPV model based on the nu-gap
% metric.  
%
% reference paper: An H?-Norm-Based Approach
% for Operating Point Selection
% and LPV Model Identification
% from Local Experiments by Dániel Vízer / Guillaume Mercère (Useful because it measures distance in an L2 sense)
% 
% Inputs : State space matrices as functions of the scheduling variable,
% minimum sceduling variable value, maximum sceduling variable value and
% threshhold for the nugap metric.
%
% Bardia Sharif
% November 2017
% Eindhoven university of technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
function [gp] = grid_nugap_metric(Af,Bf,Cf,Df,pmin,pmax,thresh)
i           =   1;
gp          =   pmin;
pstep       =   0.05*pmax;
delta_step  =   0.25*pstep;
cntp        =   0;
cntm        =   0;
while gp(i)<=pmax 
    S1          =   ss(Af(gp(i)),Bf(gp(i)),Cf(gp(i)),Df(gp(i)));
    i           =   i+1;
    gp(i)       =   gp(i-1)+ pstep;
    S2          =   ss(Af(gp(i)),Bf(gp(i)),Cf(gp(i)),Df(gp(i)));
    [~,sig]     =   gapmetric(S1,S2);    
    
    while (thresh-0.1 > sig || sig > thresh+0.1)
        if sig > thresh + 0.1
            pstep_p = pstep - delta_step;
            gp(i)   = gp(i) - pstep_p;
            cntp    = cntp+1;
        else 
            pstep_m = pstep + delta_step;
            gp(i)   = gp(i) + pstep_m;
            cntm    = cntm+1;
        end
            S2          =   ss(Af(gp(i)),Bf(gp(i)),Cf(gp(i)),Df(gp(i)));       
            [~,sig]     =   gapmetric(S1,S2);    
    end   
end
    gp(i)   =   pmax;
end