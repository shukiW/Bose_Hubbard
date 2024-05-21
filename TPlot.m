function TPlot(FockBasis,Vn)
scatter(FockBasis(:,2)-FockBasis(:,3),FockBasis(:,2)+FockBasis(:,3),30,abs(Vn).^2,'filled')
ax = gca;
ax.FontSize = 25;
colormap(jet)
axis off
end
