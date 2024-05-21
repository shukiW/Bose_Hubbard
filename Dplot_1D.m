function Dplot_1D(N,Vn,plot3D)
h1=501;
xx=(-1:2/h1:1)'; %[0,1]
h2 = 501;
yy=(0:1/h2:1)'; %[0,1]
[XX,YY]=meshgrid(xx,yy);
phi=pi*XX(:);
theta=pi*YY(:);
C = zeros(N+1,length(phi));
for m=0:N
C(m+1,:)=sqrt(nchoosek(N,m)).*(cos(theta/2)).^(N-m).*(sin(theta/2)).^m.*exp(1i*m*phi);
end
Q = abs(C'*Vn).^2;

if (plot3D)
    theta3D =reshape(theta,h2+1,h1+1);
    phi3D =reshape(phi,h2+1,h1+1);
    X = 1/2*sin(theta3D).*cos(phi3D);
    Y = 1/2*sin(theta3D).*sin(phi3D);
    Z = 1/2*cos(theta3D);
    Q3D = reshape(Q,h2+1,h1+1);
    surf(X,Y,Z,Q3D,'EdgeAlpha',0)
    view(0,0)
    xlim([-1/2 1/2])
    ylim([-1/2 1/2])
    zlim([-1/2 1/2])
else
    X = 1/2*(1-cos(theta)).*cos(phi);
    Y = 1/2*(1-cos(theta)).*sin(phi);  
    scatter(X,Y,3,Q,'filled');
    xlim([-1 1]); 
    ylim([-1 1]);
end
clim([0 0.5])
colormap(jet)
axis off square
end