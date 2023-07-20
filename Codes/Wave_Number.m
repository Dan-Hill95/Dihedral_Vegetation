function k = Wave_Number(pars,hands)
f=hands.f; g=hands.g; var=pars.sol;
syms u v mu
u0=var(1); v0=var(2); mu0=var(3);
x = -(diff(f,u) +  diff(g,v))./2;
k=sqrt(double(x(u0,v0,mu0)));
end