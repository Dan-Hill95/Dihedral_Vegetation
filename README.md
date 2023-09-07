This repository contains the MATLAB files to reproduce the data and figures from 'Predicting the emergence of localised dihedral patterns in models for dryland vegetation' by Dan J. Hill (2023).

A copy of the paper can be found here: https://arxiv.org/abs/2309.02956

The script _Quantity_Generation.m_ computes the qualitative predictors P_1, P_2, P_3 and P_4 from the paper for a two-component reaction-diffusion system. To choose one of the vegetation models considered in the paper, uncomment the relevant lines in the script.

The script _Figure_Generation.m_ uses exponential time differencing codes from

https://github.com/kleefeld80/ETDRDPIF 

to simulate the dynamics of localised vegetation patterns with some dihedral symmetry. To choose one of the vegetation models and one of the localised patterns considered in the paper, uncomment the relevant lines in the script.

If you want to adapt this code for your own two-component reaction-diffusion model, just add the function handles and parameters to _Codes/Equation.m_.
