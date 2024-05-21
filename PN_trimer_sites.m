%% calculating the vector of all Fock states

% dimension of the system
dim=(N+1)*(N+2)/2;

% the populations for each Fock vector
FockBasis=zeros(dim,3);
r1 =1+[0,cumsum(N+1:-1:1)];
r2 = N:-1:0;
for rr=1:N
    FockBasis(r1(rr):r1(rr+1)-1,2)=0:r2(rr);
    FockBasis(r1(rr):r1(rr+1)-1,1)=r2(rr):-1:0;
end
FockBasis(:,3)=N-FockBasis(:,1)-FockBasis(:,2);

%% calculating the operator matrices

% n_i = a_i^d * a_i 
n1=sparse(diag(FockBasis(:,1)));
n2=sparse(diag(FockBasis(:,2)));
n3=sparse(diag(FockBasis(:,3)));
n1s = sparse(diag(FockBasis(:,1).*(FockBasis(:,1)-1)));
n2s = sparse(diag(FockBasis(:,2).*(FockBasis(:,2)-1)));
n3s = sparse(diag(FockBasis(:,3).*(FockBasis(:,3)-1)));

cr = [1,2,3];
for rr = 1:3
    % a_i^d * a_j
    ii=cr(1); jj=cr(2); mm=cr(3);
    ind1=find(FockBasis(:,jj)>0 & FockBasis(:,ii)<N);
    ind2=zeros(length(ind1),1);
    val=zeros(length(ind1),1);

    for kk=1:length(ind1)
        ind2(kk)=find( FockBasis(:,ii)==FockBasis(ind1(kk),ii)+1 & FockBasis(:,jj)==FockBasis(ind1(kk),jj)-1 &...
                            FockBasis(:,mm)==FockBasis(ind1(kk),mm) );
        val(kk)=sqrt(FockBasis(ind1(kk),jj)*(FockBasis(ind1(kk),ii)+1));
    end

    switch ii
    case 1
        a1da2=sparse(ind2,ind1,val,dim,dim);
    case 2
        a2da3=sparse(ind2,ind1,val,dim,dim);
    case 3
        a3da1=sparse(ind2,ind1,val,dim,dim);
    end
    cr = circshift(cr,-1);
end

%% scale by N

n1=n1/N;
n2=n2/N;
n3=n3/N;
n1s=n1s/N^2;
n2s=n2s/N^2;
n3s=n3s/N^2;
a1da2=a1da2/N;
a2da3=a2da3/N;
a3da1=a3da1/N;