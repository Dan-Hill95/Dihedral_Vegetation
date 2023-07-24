%% Figure generation for localised planar patterns in vegetation models
% Dan J Hill (2022) - Saarland University

% Exponential time-stepping codes edited from:
% `A second-order exponential time differencing scheme for non-linear reaction-diffusion systems with dimensional splitting'
% Asante-Asamani et al. (2020) - https://github.com/kleefeld80/ETDRDPIF

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

%% Choose initial guess for D_m pattern

%%% Hexagon %%%

amp.x = (1/3)*[1,1,1];
amp.m = 6;
if ProbClass=='Kl'
amp.eps = 0.00005;
elseif ProbClass=='KG'
amp.eps = 0.0005;
elseif ProbClass=='vH'
if ProbName=='vH_1'
amp.eps = 0.0002;
    elseif ProbName=='vH_2'
amp.eps = 0.0002;
    end
elseif ProbClass=='Gi'
amp.eps = 0.0005;
end
SolnClass= 'Hex';


%%% Square %%%

% amp.x = (1/6)*[-2,1,1,-2,1,1];
% amp.m = 4;
% if ProbClass=='Kl'
% amp.eps = 0.00002;
% elseif ProbClass=='KG'
% amp.eps = 0.0002;
% elseif ProbClass=='vH'
% if ProbName=='vH_1'
% amp.eps = 0.0002;
%     elseif ProbName=='vH_2'
% amp.eps = 0.0002;
%     end
% elseif ProbClass=='Gi'
% amp.eps = 0.0002;
% end
% SolnClass= 'Sqr';

%%% Pentagon %%%

% amp.x = (1/3)*[-1,1,1,1];
% amp.m = 5;
% amp.eps = 0.0001;
% SolnClass= 'Pnt';
% amp.x = (1/3)*[-1,1,1,1];
% amp.m = 5;
% if ProbClass=='Kl'
% amp.eps = 0.00002;
% elseif ProbClass=='KG'
% amp.eps = 0.0002;
% elseif ProbClass=='vH'
% if ProbName=='vH_1'
% amp.eps = 0.0002;
%     elseif ProbName=='vH_2'
% amp.eps = 0.0002;
%     end
% elseif ProbClass=='Gi'
% amp.eps = 0.0002;
% end
% SolnClass= 'Pnt';



%%% Triangle %%%

% amp.x = (1/6)*[-1,sqrt(3),1,-sqrt(3),1,-sqrt(3)];
% amp.m = 3;
% if ProbClass=='Kl'
% amp.eps = 0.001;
% elseif ProbClass=='KG'
% amp.eps = 0.001;
% elseif ProbClass=='vH'
%     if ProbName=='vH_1'
% amp.eps = 0.0001;
%     elseif ProbName=='vH_2'
% amp.eps = 0.002;
%     end
% elseif ProbClass=='Gi'
% amp.eps = 0.001;
% end
% SolnClass= 'Tri';



%%% Super-Hexagon %%%

% amp.x = (1/3)*[1,1,-1,-1];
% amp.m = 6;
% if ProbClass=='Kl'
% amp.eps = 0.0000002;
% elseif ProbClass=='KG'
% amp.eps = 0.0002;
% elseif ProbClass=='vH'
% if ProbName=='vH_1'
% amp.eps = 0.0002;
%     elseif ProbName=='vH_2'
% amp.eps = 0.0007;
% end
% elseif ProbClass=='Gi'
% amp.eps = 0.0002;
% end
% SolnClass= 'SHx';



%%% Dodecagon %%%

% amp.x = (1/3)*[1,1,1];
% amp.m = 12;
% amp.eps = 0.000001;
% SolnClass= 'Ddc';



%%% Heptagon %%%

% amp.x = (1/3)*[-1,1,1,1];
% amp.m = 7;
% amp.eps = 0.0001;
% SolnClass= 'Hpt';


%%% Nonagon %%%

% amp.x = (1/3)*[-1,1,1];
% amp.m = 9;
% amp.eps = 0.0001;
% SolnClass= 'Non';


%%% Radial Spot %%%

% amp.x = (1/1)*[1,0,0];
% amp.m = 0;
% amp.eps = 0.0001;
% SolnClass= 'Rad';


%%% Rhombus %%%

% amp.x = (10)*[1,1,1];
% amp.m = 2;
% amp.eps = 0.0001;
% SolnClass= 'Rmb';


%%% Alternate Square %%%

% amp.x = (10)*[1,1,1];
% amp.m = 4;
% % amp.eps = 0.0001;
% % amp.eps = 0.001;
% amp.eps = 0.005;
% SolnClass= 'Ssq';

%% Initial Guess and Time Stepping

% Name solution folder
FolderName = [ProbName,'_',SolnClass];
if isfolder(FolderName)==1
rmdir(FolderName,'s');
end
mkdir(FolderName)

% Solve Galerkin matching equation for D_m pattern and define predicted localised solution
[Sol.x, Sol.y, Sol.usol, Sol.vsol, Sol.musol] = Initial_Guess(amp,pars);

if ProbClass == 'vH'
if ProbName == 'vH_1'
temp1=Sol.usol - Sol.usol(end);
temp2=Sol.vsol - Sol.vsol(end);
Sol.usol = 30*temp1 + Sol.usol(end);
Sol.vsol = 30*temp2 + Sol.vsol(end);
end
end

% Implement exponential time-stepper for initial guess
[ts.runtime,ts.w_old]=Time_Stepper(pars,Sol,hands,FolderName);

% Figures saved in solution folder at every t= n*100, video saved after completion
close all

% Clear extra structures
clear("FolderName","ProbName","ProbClass","SolnClass","Init","options","prob","ts")
