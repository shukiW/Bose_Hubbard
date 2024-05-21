function plot_dimer_seperatrix(N,U,J,Delta,plot3D)

u = U*N;
pol = [u^2, -Delta*u, (J^2+Delta^2-u^2)/4, Delta*u/4, -Delta^2/16];
r = roots(pol);
r = r(abs(imag(r))<10^-5);
CSx = -J*r./(2*u*r-Delta);
CSz =  r;

E_s = sort(-J*CSx -Delta*CSz +u*CSz.^2);
Ej = E_s(2);
nn = linspace(0,1,10^5);
Csphi = (u*(1/2-nn).^2-Delta*(1/2-nn)-Ej)./(J*sqrt((1-nn).*nn));
msk3 = Csphi <= 1 & Csphi >= -1;
Phi1 = acos(Csphi);
Phi2 = -acos(Csphi);

hold on
    if (plot3D)
        scatter3(sqrt(nn(msk3).*(1-nn(msk3))).*cos(Phi1(msk3)),sqrt(nn(msk3).*(1-nn(msk3))).*sin(Phi1(msk3)),1/2-nn(msk3),2,'w','filled')
        scatter3(sqrt(nn(msk3).*(1-nn(msk3))).*cos(Phi2(msk3)),sqrt(nn(msk3).*(1-nn(msk3))).*sin(Phi2(msk3)),1/2-nn(msk3),2,'w','filled')
    else
        scatter(nn(msk3).*cos(Phi1(msk3)),nn(msk3).*sin(Phi1(msk3)),2,'w','filled')
        scatter(nn(msk3).*cos(Phi2(msk3)),nn(msk3).*sin(Phi2(msk3)),2,'w','filled')
    end

end