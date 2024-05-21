%% calculate the vector of all Fock states

% dimension of the system
dim=N+1;
% the populations for each Fock vector
FockBasis=zeros(dim,2);

FockBasis(:,2)=0:N;
FockBasis(:,1)=N-FockBasis(:,2);

%% calculate the operator matrices

% n_i = a_i^d * a_i 
n1=sparse(diag(FockBasis(:,1)));
n2=sparse(diag(FockBasis(:,2)));
n1s = sparse(diag(FockBasis(:,1).*(FockBasis(:,1)-1)));
n2s = sparse(diag(FockBasis(:,2).*(FockBasis(:,2)-1)));

% a_i^d * a_j
ii=1; jj=2;
ind1=find(FockBasis(:,jj)>0 & FockBasis(:,ii)<N);
ind2=zeros(length(ind1),1);
val=zeros(length(ind1),1);

for kk=1:length(ind1)
    ind2(kk)=find( FockBasis(:,ii)==FockBasis(ind1(kk),ii)+1 & FockBasis(:,jj)==FockBasis(ind1(kk),jj)-1);
    val(kk)=sqrt(FockBasis(ind1(kk),jj)*(FockBasis(ind1(kk),ii)+1));
end

a1da2=sparse(ind2,ind1,val,dim,dim);

% a_i^d * a_i^d * a_j * a_j
ii=1; jj=2;
ind1=find(FockBasis(:,jj)>1 & FockBasis(:,ii)<(N-1));
ind2=zeros(length(ind1),1);
val=zeros(length(ind1),1);

for kk=1:length(ind1)
    ind2(kk)=find( FockBasis(:,ii)==FockBasis(ind1(kk),ii)+2 & FockBasis(:,jj)==FockBasis(ind1(kk),jj)-2);
    val(kk)=sqrt(FockBasis(ind1(kk),jj)*(FockBasis(ind1(kk),jj)-1)*(FockBasis(ind1(kk),ii)+1)*(FockBasis(ind1(kk),ii)+2));
end

a1da1da2a2=sparse(ind2,ind1,val,dim,dim);


%------------------------------------------------------------
n1=n1/N;
n2=n2/N;
n1s=n1s/N^2;
n2s=n2s/N^2;
a1da2=a1da2/N;
a1da1da2a2 = a1da1da2a2/N^2;
%-------------------------------------------------------------