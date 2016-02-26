%plot displacement
function plotfun(COORD,vm,data,nameplot,num)
figure(num)
set(num,'Position',[20+400*(num-1) 20 400 400]);
vmM=max(vm);
xyplot=linspace(-0.03,0.03,200);
[Xp,Yp]=meshgrid(xyplot,xyplot);
[Xr,Yr,Fr] = griddata(COORD(:,1),COORD(:,2),vm,Xp,Yp);
%vlevel=[0:.2:2];
%contourf(Xr,Yr,Fr,vlevel)
contour(Xr,Yr,Fr,20)
%colorbar
if isempty(data.hole_nodes1)
    hole=[];
else
    hole=convhull(COORD(data.hole_nodes1,:));
end
hold on
patch(COORD(data.hole_nodes1(hole),1),COORD(data.hole_nodes1(hole),2),vmM*(hole./hole)','white')
if isempty(data.hole_nodes2)
    hole=[];
else
    hole=convhull(COORD(data.hole_nodes2,:));
end
hold on
patch(COORD(data.hole_nodes2(hole),1),COORD(data.hole_nodes2(hole),2),vmM*(hole./hole)','white')
set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
set(gca,'Visible','off')
%ti=['v=' num2str(vfrac) ';   g_{vM} = ' num2str(theta) ';   it = ' num2str(it)]
%title(ti)
axis equal
axis tight
set(gcf,'paperunits','centimeters')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'papersize',[8,8])
set(gcf,'paperposition',[0,0,8,8])
saveas(gcf,[nameplot num2str(num) '.eps'],'psc2');
end

