function [F,J] = Turing_Conditions(var,hands)
%% Function F
f=hands.f;g=hands.g;
%Computing as symbolic functions
syms u v mu
u0=var(1); v0=var(2); mu0=var(3);
x = 4*diff(f,v).*diff(g,u) + (diff(f,u) -  diff(g,v)).^2;

%Function - output as double
F = [double(f(u0,v0,mu0)); double(g(u0,v0,mu0)); double(x(u0,v0,mu0))];
% F = [f ; g ; x]

%% Jacobian J
if nargout>1

%Partial derivatives - Comuting as symb  
fu=diff(f,u);
fv=diff(f,v);
fmu=diff(f,mu);
gu=diff(g,u);
gv=diff(g,v);
gmu=diff(g,mu);
xu=diff(x,u);
xv=diff(x,v);
xmu=diff(x,mu);

%Jacobian - output as double
J = [double(fu(u0,v0,mu0)),double(fv(u0,v0,mu0)),double(fmu(u0,v0,mu0));...
    double(gu(u0,v0,mu0)),double(gv(u0,v0,mu0)),double(gmu(u0,v0,mu0));...
    double(xu(u0,v0,mu0)),double(xv(u0,v0,mu0)),double(xmu(u0,v0,mu0))];
end
end