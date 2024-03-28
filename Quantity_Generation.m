%% Quantity generation for localised planar patterns in vegetation models
% Dan J Hill (2022) - Saarland University
function [quant,pars]=Quantity_Generation(ProbClass,Init)
close all 
clc

addpath Codes
% [In order to check individual matlab functions, see the "Codes" folder]

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
quant.P4=pars.c3;

% Clear extra structures
clear("FolderName","ProbClass","SolnClass","Init","options","prob","ts","hands")

%% Print Results 

% Print each predictor
fprintf('P1=%d\n',quant.P1);
fprintf('P2=%d\n',quant.P2);
fprintf('P3=%d\n',quant.P3);
fprintf('P4=%d\n',quant.P4);

% Print interpretation of results

% Bifurcation direction
if quant.P1<0
fprintf('Localised patterns emerge for %c < %d\n',956,pars.sol(3));
elseif quant.P1>0
fprintf('Localised patterns emerge for %c > %d\n',956, pars.sol(3));
else
    fprintf('Inconclusive: P1 is zero\n');
end

% Profile of spot A-type patterns
if quant.P2<0
phase = 'in-phase';
elseif quant.P2>0
phase = 'anti-phase';
end
if quant.P3<0
polarity = 'gaps';
elseif quant.P3>0
polarity = 'peaks';
end

if quant.P2==0
fprintf('Inconclusive: P2 is zero\n');
elseif quant.P3==0
fprintf('Inconclusive: P3 is zero\n');
else
fprintf(['Spot patterns emerge as ',phase,' ',polarity,'\n']);
end

% Emergence of ring patterns
if quant.P4<0
fprintf('Rings and unstable stripes do emerge\n');
elseif quant.P4>0
fprintf('Rings and unstable stripes do not emerge\n');
else
    fprintf('Inconclusive: P4 is zero\n');
end
end