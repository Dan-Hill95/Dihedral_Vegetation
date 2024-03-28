function [pars,hands] = Equation(ProbClass)
syms u v mu

% pars structure contains parameters
% hands structure contains function handles


if ProbClass == 'My' % User defined model

    if exist('My_Model.mat')==0
        error('No model defined. Run my_model.m first.');
    else
load('My_Model.mat',"model");
    end
p = model.p;
q = model.q;
f0 = model.f0;
g0 = model.g0;

elseif ProbClass == 'Kl' % Logistic Klausmeier Model

p = [0.45,1];
q=[182.5,0];

%Parameters
m=p(1);b=p(2);

f0(u,v,mu) = -(1-b*u).*v.*u.^2 + m.*u;
g0(u,v,mu) = -(mu-v-v.*u.^2);

elseif ProbClass == 'KG' % Klausmeier-Gray-Scott Model

p = [0.5];
q=[7.2,0];

%Parameters
m=p(1);

f0(u,v,mu) = -v.*u.^2 + m.*u;
g0(u,v,mu) = -(mu-v-v.*u.^2);


elseif ProbClass == 'vH' % von Hardenberg Model 


p = [1.6, 1.6, 0.2, 1.5];
q=[100,3];

%Parameters
ga=p(1);s=p(2);n=p(3);r=p(4);


f0(u,v,mu) =-((ga.*v)./(1+s.*v)-n-u).*u;
g0(u,v,mu) =-(mu-(1-r.*u).*v - u.*v.^2);

elseif ProbClass == 'Gi' % NFC-Gilad Model

%Zelnik et al. (2015) paper    
p = [16/35, 14/5, 10/7, 7/10];
q=[125,0];
    
%Parameters
L=p(1);e=p(2);n=p(3);r=p(4);


f0(u,v,mu) = -L.*v.*u.*(1-u).*(1+e.*u).^2 + u;
g0(u,v,mu) = -(mu - n.*(1-r.*u).*v - L.*u.*v.*(1+e.*u).^2);

else
    error('Problem class not recognised.')
end

Dv = q(1); beta = q(2);

%Function RHS's for Steady States
f = f0;
g = (1/Dv).*g0 + beta*f0;         

%Function RHS's for time-stepping
if Dv ==1
fx = matlabFunction(f0(u,v, mu),'Vars',[mu,u,v]);
gx = matlabFunction(g0(u,v, mu),'Vars',[mu,u,v]);
else
fx = matlabFunction(f0(u,v+(beta*Dv/(Dv-1)).*u, mu),'Vars',[mu,u,v]);
gx = matlabFunction(g0(u,v+(beta*Dv/(Dv-1)).*u, mu) - (beta*Dv/(Dv-1)).*f0(u,v+(beta*Dv/(Dv-1)).*u, mu),'Vars',[mu,u,v]);
end
pars.p=p;pars.q=q;
hands.f=f; hands.g=g;
hands.f0=fx;hands.g0=gx;

end