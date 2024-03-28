% Define your own two-commponent reaction-diffusion model of the form
% 
%           u_t =  \Delta u - f0(u,v,mu)
%           v_t = D\Delta (v - \beta u) - g0(u,v,mu)

syms u v mu

p = []; % Fixed reaction parameters                   e.g. p = [0.45,1];
q = [ , ]; % Diffusion parameters q = [delta, beta]     e.g. q = [182.5,0];

% Now define the reaction functions 

f0(u,v,mu) = ;          % e.g. f0(u,v,mu) = -(1-p(2).*u).*v.*u.^2 + p(1).*u;
g0(u,v,mu) = ;          % e.g. g0(u,v,mu) = -(mu-v-v.*u.^2);


% (Make sure to use element-wise multiplication .*, 
% so that f0 and g0 can be converted to vector functions later)

% Save parameters and functions in one structure 

model.p = p; model.q =q; model.f0 =f0; model.g0 = g0;
save("My_Model.mat",'model','-mat');

clear all