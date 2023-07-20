%Dan J Hill (2021) - Solving the simple D_{2k} patch matching condition
%x is the initial guess: an (N+1)-dimensional vector of the form x = [x(1), x(2), ..., x(N+1)]
%k is the lattice index: solutions are invariant under rotations of pi/k
%r_max is the maximum range of the mesh
%mu is the localisation of the solution

%For 3|k, a good initial guess will be of the form x=y*[1,...,1], for some
%small y

%For 3~|k, a good initial guess will be of the form x=y*[-0.5, 1, 1, -0.5,1,1,...,-0.5,1,1], for some
%small y
function a_out=MatchSoln(x,k)
%Introduce parameters
N=length(x)-1;
a0=x.*ones(1,N+1);
%fsolve options
options = optimset('Display','iter','TolFun',1e-9);
%Solving the (N+1) D_{k} matching condition
a_out=fsolve(@(a) match(a,k),a0,options);

end