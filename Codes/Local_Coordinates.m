function [M1,M2,Q,C,U0,V0,U1,V1] = Local_Coordinates(pars,hands)
f=hands.f;g=hands.g;k=pars.k;var=pars.sol;
syms u v mu  u1 v1 u2 v2 %Symbolic variables - for computations
u0=var(1); v0=var(2); mu0=var(3);   %Double variables   - for output

%% Derivative of Uniform state
S=solve([f==0; g==0],[u,v],'MaxDegree',3);
vs(mu)=S.v;
vss=double(vs(mu0));
for i= 1:length(vss)
if abs(vss(i) - v0)<1e-06
 ux = S.u(i);
 vx = S.v(i);
end
end
%% Symbolic function derivatives

fu      =   diff(f,u);
fv      =   diff(f,v);
fuu     =   diff(f,u,u);
fvu     =   diff(f,u,v); 
fux     =   fu(ux,vx,mu);
fvx     =   fv(ux,vx,mu);
fmuu(u,v,mu)=diff(fux,mu);
fmuv(u,v,mu)=diff(fvx,mu);
fvv     =   diff(f,v,v);
fuuu    =   diff(f,u,u,u);
fuuv    =   diff(f,v,u,u);
fuvv    =   diff(f,v,v,u);
fvvv    =   diff(f,v,v,v);

gu      =   diff(g,u);
gv      =   diff(g,v);
guu     =   diff(g,u,u);
gvu     =   diff(g,u,v);
gux     =   gu(ux,vx,mu);
gvx     =   gv(ux,vx,mu);
gmuu(u,v,mu)=diff(gux,mu);
gmuv(u,v,mu)=diff(gvx,mu);
gvv     =   diff(g,v,v);
guuu    =   diff(g,u,u,u);
guuv    =   diff(g,v,u,u);
guvv    =   diff(g,v,v,u);
gvvv    =   diff(g,v,v,v);

%% Eigenfunctions

U0 = [double(fv(u0,v0,mu0)); -(k^2 + double(fu(u0,v0,mu0)))];
V0 = [1./(double(fv(u0,v0,mu0)));0];
U1 = [0 ; k^2];
V1 = (1./((k^2).*double(fv(u0,v0,mu0))))*[(k^2 + double(fu(u0,v0,mu0))) ; double(fv(u0,v0,mu0))];

%% Linear Terms
M1 = [double(fu(u0,v0,mu0)), double(fv(u0,v0,mu0)); double(gu(u0,v0,mu0)), double(gv(u0,v0,mu0))];
M2 = real([double(fmuu(u0,v0,mu0)), double(fmuv(u0,v0,mu0)); double(gmuu(u0,v0,mu0)), double(gmuv(u0,v0,mu0))]);

%% Quadratic Terms
%Quadratic function
F2(u,v,mu) = (1/2)*[double(fuu(u0,v0,mu0)), double(fvu(u0,v0,mu0)), double(fvv(u0,v0,mu0)); ...
    double(guu(u0,v0,mu0)), double(gvu(u0,v0,mu0)), double(gvv(u0,v0,mu0))]*[u.^2 ; 2.*u.*v; v.^2];
%Convert into bilinear function
Q(u,v,u1,v1) = (1/2)*(F2(u+u1,v+v1,mu0) - F2(u,v,mu0) - F2(u1,v1,mu0));

%% Cubic Terms
%Cubic function
F3(u,v,mu) = (1/6)*[double(fuuu(u0,v0,mu0)), double(fuuv(u0,v0,mu0)), double(fuvv(u0,v0,mu0)), double(fvvv(u0,v0,mu0)); ...
    double(guuu(u0,v0,mu0)), double(guuv(u0,v0,mu0)), double(guvv(u0,v0,mu0)), double(gvvv(u0,v0,mu0))...
    ]*[u.^3 ; 3.*u.^2.*v; 3.*u.*v.^2; v.^3];
%Convert into trilinear function
C(u,v,u1,v1,u2,v2) = (1/6)*(F3(u+u1+u2,v+v1+v2,mu0)+ F3(u,v,mu0) + F3(u1,v1,mu0) +F3(u2,v2,mu0) ...
    - F3(u+u1,v+v1,mu0) - F3(u+u2,v+v2,mu0) - F3(u1+u2,v1+v2,mu0));

end