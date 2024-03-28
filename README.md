This repository contains the MATLAB files to reproduce the data and figures from 'Predicting the emergence of localised dihedral patterns in models for dryland vegetation' by Dan J. Hill (2023).

A copy of the paper can be found here: https://arxiv.org/abs/2309.02956

The script _Figure_Generation.m_ uses exponential time differencing codes from

https://github.com/kleefeld80/ETDRDPIF 

to simulate the dynamics of localised vegetation patterns with some dihedral symmetry. To choose one of the vegetation models and one of the localised patterns considered in the paper, uncomment the relevant lines in the script.

If you want to adapt this code for your own two-component reaction-diffusion model, just add the function handles and parameters to _Codes/Equation.m_.

--------------------------------------------------------------------------------------------------------------------------------------------------------------

The script _Quantity_Generation.m_ computes the qualitative predictors P_1, P_2, P_3 and P_4 from the paper for a two-component reaction-diffusion system. To run this script, you write

Quantity_Generation(ProbClass,Init);

where *ProbClass* is a two-letter string that defines the reaction-diffusion model and *Init* is an initial guess of a Turing point \[u*, v*, &#956*\].

For your choice of model, you can either define your own reaction-diffusion model in _my_model.m_ and select ProbClass = 'My' with an appropriate initial guess Init, or choose one of the pre-defined vegetation models listed below:

**Klausmeier-Gray-Scott model:**
ProbClass = 'KG', Init = \[1.0708, 0.4669 , 1.0023\];

**Logistic Klausmeier model:**
ProbClass = 'Kl', Init = \[0.465, 1.809,2.200\];

**Gilad et al. model:** 
ProbClass = 'Gi', Init = \[0.4743,0.7678,1.6350\];

**von Hardenberg model (Turing point 1):** 
ProbClass = 'vH', Init = \[0.0165, 0.173, 0.169\];

**von Hardenberg model (Turing point 2):** 
ProbClass = 'vH', Init = \[0.271, 0.556, 0.414\];
