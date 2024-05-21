clear;
tStart=tic;
N=10;
PN_dimer;
Sz = N/2*(n1-n2);
S2z = Sz^2;
Sx = N/2*(a1da2+(a1da2)');
toc(tStart)

%%
U0 = 1/N;
J0 = 2;
J1 = -2;
Delta0 = 0;
VJ = 1/10;
h=3000;
J = linspace(J1,J0,h+1);
t = J/VJ;
s = length(t);

E = zeros(dim,s);
Prob = zeros(dim,s);
for k=1:s
    % the Hamiltonian
    H = -Delta0*Sz -J(k)*Sx +U0*S2z;
    [V,D]=eig(full(H));
    if (k==1)
        % initial state
        Vn = V(:,1);
    else
        % time evolution
        Vn = expm(-1i*full(H)*(t(k)-t(k-1)))*Vn;
    end        
    E(:,k)=diag(D);
    Prob(:,k) = abs(V'*Vn).^2;  
end
fprintf('%d\n',k);
toc(tStart)

% probability distribution
figure()
ttable = ones(dim,1)*t;
msk = Prob > 0.04;
scatter(ttable(~msk),E(~msk),2,'MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5]);
hold on
scatter(ttable(msk),E(msk),6,Prob(msk),'filled');
colormap(jet)
xlabel('$t$','Interpreter','latex');
ylabel('$E$','Interpreter','latex');
clim([0 1])
xlim([t(1) t(end)])
ax = gca;
ax.FontSize = 25;