function [x, y, usol, vsol, musol] = Initial_Guess(amp,pars)
a=amp.x;
m=amp.m;
k=pars.k;
var=pars.sol;
eps0=amp.eps;
U0=pars.U0;
c0=pars.c0;
gamma=pars.gamma;

u0=var(1);v0=var(2);mu0=var(3);

Ln = ceil(20*pi/k);

a0 = MatchSoln(a,m);
n=length(a0);
N = Ln;
r = 0:sqrt(2)*Ln/(N-1):sqrt(2)*Ln;

for i=0:n-1
u1(1+(i)*N:(i+1)*N)= ((2*k*sqrt(3))./gamma).*u0.*a0(i+1).*besselj(m*(i),k*r).*U0(1).*exp(-sqrt(eps0).*r);
v1(1+(i)*N:(i+1)*N)= ((2*k*sqrt(3))./gamma).*v0.*a0(i+1).*besselj(m*(i),k*r).*U0(2).*exp(-sqrt(eps0).*r);
end

t = 0:0.01:2*pi;
%    t = t';
  u = u1(1:n*N);
  v = v1(1:n*N);
  [R,T]=meshgrid(r,t);
  [U,~]=meshgrid(u,t);
  [V,~]=meshgrid(v,t);
  UU(:,:,1)=U(:,1:N)/2;
  VV(:,:,1)=V(:,1:N)/2;
  for j=1:n-1
      UU(:,:,j+1)= U(:,1+j*N:(j+1)*N);
      UU(:,:,j+1)= UU(:,:,j+1).*cos(m*j.*T);
      VV(:,:,j+1)= V(:,1+j*N:(j+1)*N);
      VV(:,:,j+1)= VV(:,:,j+1).*cos(m*j.*T);
  end
  U0 = sum(UU,3);
  V0 = sum(VV,3);
X =R.*cos(T);
Y =R.*sin(T);
xRange = linspace(-Ln,Ln,N);
yRange = linspace(-Ln,Ln,N);
[x,y] = meshgrid(xRange,yRange);
warning('off','all');
ux= griddata(X,Y,U0,x,y);
vx= griddata(X,Y,V0,x,y);
warning('on','all');
ux=fillmissing(ux,'constant',0);
vx=fillmissing(vx,'constant',0);

usol=u0+ux;
usol = reshape(usol,[],1);
vsol=v0+vx;
vsol = reshape(vsol,[],1);
musol = mu0 + eps0./c0;
end