%% Quantity generation for localised planar patterns in vegetation models
% Dan J Hill (2022) - Saarland University

close all 
clear
clc

addpath Codes
% [In order to check individual matlab functions, see the "Codes" folder]

%% Choose model - Problem class & initial guess for Turing points

%%% Klausmeier %%%

% ProbClass='Kl';
% Init=[0.465, 1.809,2.200];
% ProbName='Kl';


%%% Klausmeier-Gray-Scott %%%

% ProbClass='KG';
% Init=[ 1.0708, 0.4669 , 1.0023];
% ProbName='KG';


%%% von Hardenberg - Turing point 1 %%%

% ProbClass='vH';
% Init = [0.0165, 0.173, 0.169];
% ProbName='vH_1';


%%% von Hardenberg - Turing point 2 %%%

ProbClass='vH';
Init = [0.271, 0.556, 0.414];
ProbName='vH_2';

%%% NFC - Gilad %%%

% ProbClass='Gi';
% Init=[0.4743,0.7678,1.6350];
% ProbName='Gi';

%% Equations, Turing points, conditions

% Input problem class to generate parameters and function handles (for both steady states and time-stepping)
[pars,hands]=Equation(ProbClass);

% Solve algebraic conditions for a Turing bifurcation of a uniform state, near some initial guess
prob = @(var)Turing_Conditions(var,hands);
options = optimset('Jacobian','on','Display','iter','MaxIter',50,'TolFun',1e-7,'DerivativeCheck','off');
pars.sol = fsolve(prob,Init,options);

% Compute wave number k
pars.k = Wave_Number(pars,hands);

% Check k is real
if abs(imag(pars.k)) > 1e-04
    error('Turing point not found - wave number is complex')
end

% Compute local operators for spatial dynamics
[hands.M1,hands.M2,hands.Q,hands.C,pars.U0,pars.V0,pars.U1,pars.V1] = Local_Coordinates(pars,hands);

% Compute constants c0, gamma & c3
[pars.c0,pars.gamma,pars.c3] = Localised_Conditions(pars,hands);

% Construct structure of qualitative predictors
quant.P1=pars.c0;
quant.P2=pars.U0(2)./pars.U0(1);
quant.P3=pars.U0(1)./pars.gamma;
quant.P4=pars.c3/pars.c0;

% Clear extra structures
clear("FolderName","ProbName","ProbClass","SolnClass","Init","options","prob","ts","hands")

% Print each predictor
fprintf('P_1=%d\n',quant.P1);
fprintf('P_2=%d\n',quant.P2);
fprintf('P_3=%d\n',quant.P3);
fprintf('P_4=%d\n',quant.P4);