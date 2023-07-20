function [pars,hands] = Equation(ProbClass)
syms u v mu

% pars structure contains parameters
% hands structure contains function handles


if ProbClass == 'Kl' % Klausmeier Model

p = [0.45,1];
q=[182.5,0];

%Parameters
m=p(1);b=p(2);
%Approximate Turing Point: [ 0.4667, 2.009, 2.446]

f0(u,v,mu) = -(1-b*u)*v.*u.^2 + m.*u;
g0(u,v,mu) = -(mu-v-v.*u.^2);

elseif ProbClass == 'KG' % Klausmeier-Gray-Scott Model

p = [0.5];
q=[7.2,0];

%Parameters
m=p(1);
%Approximate Turing Point: [ 0.4667, 2.009, 2.446]

f0(u,v,mu) = -v.*u.^2 + m.*u;
g0(u,v,mu) = -(mu-v-v.*u.^2);


elseif ProbClass == 'vH' % von Hardenberg Model 


p = [1.6, 1.6, 0.2, 1.5];
q=[100,3];

%Parameters
ga=p(1);s=p(2);n=p(3);r=p(4);
%Approximate Turing Point: [0.0165, 0.173, 0.169]
%Approximate Turing Point: [0.271, 0.556, 0.414]


f0(u,v,mu) =-((ga.*v)./(1+s.*v)-n-u).*u;
g0(u,v,mu) =-(mu-(1-r.*u).*v - u.*v.^2);


elseif ProbClass == 'Br' % Brusselator Model

p = [16.421];
q=[1/(0.01)^2,0];

%Parameters
a=p(1);
%Approximate Turing Point: [16.421, 0.0825, 1.3554]

f0(u,v,mu) = -(a- (mu+1).*u + v.*u.^2);
g0(u,v,mu) = -mu*u + v.*u.^2;


elseif ProbClass == 'Gi' % NFC-Gilad Model

%Zelnik et al. (2015) paper    
p = [16/35, 14/5, 10/7, 7/10];
q=[125,0];
%Approximate Turing Point: [0.4743,0.7678,1.6350]

%Zelnik et al. (2016) paper    
% p = [2, 6, 2, 0.2];
% q=[1000,0];
%Approximate Turing Point: [0.5265,0.0611,1.2210]
    
%Parameters
L=p(1);e=p(2);n=p(3);r=p(4);


f0(u,v,mu) = -L.*v.*u.*(1-u).*(1+e.*u).^2 + u;
g0(u,v,mu) = -(mu - n.*(1-r.*u).*v - L.*u.*v.*(1+e.*u).^2);

elseif ProbClass == 'SH' % Swift-Hohenberg Model

p = [1.6,-1];
q=[1,0];
    
%Parameters
ga=p(1);ka=p(2);
%Approximate Turing Point: [0,0,0]

f0(u,v,mu) = -u + v;
g0(u,v,mu) = -v - mu.*u + ga.*u.^2 + ka.*u.^3;


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