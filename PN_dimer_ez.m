%% calculate the vector of all Fock states

% dimension of the system
dim=N+1;
% the populations for each Fock vector
FockBasis=zeros(dim,2);

FockBasis(:,2)=0:N;
FockBasis(:,1)=N-FockBasis(:,2);

%% calculate the operator matrices
j = N/2; m = -j:j-1;
Jp = sparse(diag(sqrt(j*(j+1)-m.*(m+1)),1));
Jm = Jp';
Jx = (Jm+Jp)/2;
Jy = 1i*(Jm-Jp)/2;

n0 = sparse(diag(FockBasis(:,1)));
n1 = sparse(diag(FockBasis(:,2)));
Jz = 1/2*(n0-n1);
% note that Jz is also Jz = (Jp*Jm-Jm*Jp)/2; 
