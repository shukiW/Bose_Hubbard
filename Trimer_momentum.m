clear;
tStart=tic;

N = 60;
u = 2.3;
flux = 0.2*pi;
Delta = 10^-6;
PN_trimer_momentum;
Bias = 1i*Delta/sqrt(3)*(b1db2+b2db3+b3db1-b1db2'-b2db3'-b3db1');
Ints = n1s+n2s+n3s;
Intp = 4*n1*n2 + 4*n2*n3 + 4*n3*n1; 
IntH0 = Ints + Intp + 2*b1db1db2b3 + 2*b1db1db2b3';
IntHmp = 2*b2db2db3b1 + 2*b2db2db3b1' + 2*b3db3db1b2 + 2*b3db3db1b2';
toc(tStart)

E1 = -cos(flux/3);
E2 = -cos(2*pi/3 - flux/3);
E3 = -cos(4*pi/3 - flux/3);
H = Bias + E1*n1+E2*n2+E3*n3 + u/6*(IntH0+IntHmp);
[V,D]=eig(full(H));

Vi = V(:,end-5);

figure()
TPlot(FockBasis,Vi)