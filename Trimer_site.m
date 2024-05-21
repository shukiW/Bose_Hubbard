clear;
tStart=tic;

N = 60;
u = 2.3;
flux = 0.2*pi;
Delta = 10^-6;
PN_trimer_sites;
Bias =  Delta*n2 - Delta*n3; 
W = a1da2*exp(-1i*flux/3)+(a1da2')*exp(1i*flux/3)...
      +a2da3*exp(-1i*flux/3)+(a2da3')*exp(1i*flux/3)...
      +a3da1*exp(-1i*flux/3)+(a3da1')*exp(1i*flux/3);
Int = n1s+n2s+n3s;
toc(tStart)

H = Bias + u/2*Int - 1/2*W;
[V,D]=eig(full(H));
   
Vi = V(:,end-5);

figure()
TPlot(FockBasis,Vi)