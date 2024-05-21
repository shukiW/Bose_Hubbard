clear;
tStart=tic;

N = 30;
PN_dimer_ez;
toc(tStart)

u = 5;
U = 1/N;
J = U*N./u;
Delta = 10^-6*U;


H = -J*Jx -Delta*Jz +U*Jz^2;
[V,D] = eig(full(H));

Vi = V(:,17);

figure()
Dplot_1D(N,Vi,0);
plot_dimer_seperatrix(N,U,J,Delta,0)

figure()
Dplot_1D(N,Vi,1);
plot_dimer_seperatrix(N,U,J,Delta,1)