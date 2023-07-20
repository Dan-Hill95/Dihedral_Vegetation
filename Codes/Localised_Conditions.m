function [c0,gamma,c3] = Localised_Conditions(pars,hands)
M2=hands.M2;
Q=hands.Q;
C=hands.C;
U0=pars.U0;
V0=pars.V0;
U1=pars.U1;
V1=pars.V1;

Q00 = double(Q(U0(1),U0(2),U0(1),U0(2)));
Q01 = double(Q(U0(1),U0(2),U1(1),U1(2)));
C000 = double(C(U0(1),U0(2),U0(1),U0(2),U0(1),U0(2)));

c0 = (1/4)*dot(V1, -M2*U0);
gamma = dot(V1,Q00);
c3 = -((5/6)*(dot(V0,Q00)+ dot(V1,Q01)) + (19/18)*dot(V1,Q00))*dot(V1,Q00) - (3/4)*dot(V1,C000);
end