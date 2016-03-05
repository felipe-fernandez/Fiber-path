%plot scalar nodal values fieldp
function plotfun(COORD,fieldp,dene,data,nameplot,num,lay,vel)
figure(num)
set(num,'Position',[20+400*(num-1) 20 400 400]);
vmM=max(fieldp);
xyplot=linspace(-0.03,0.03,200);
[Xp,Yp]=meshgrid(xyplot,xyplot);

if (num~=(data.nlay+1) || lay==1) && vel~='c'
    [Xpd,Ypd,Fd] = griddata(data.xc_el(:,1),data.xc_el(:,2),dene,Xp,Yp);
    colormap(flipud(gray))
    contourf(Xpd,Ypd,Fd,[0 .25 .5 .75 1],'edgecolor','none');
    caxis([0 1])
    %surf(Xpd,Ypd,Fd,'FaceAlpha',0.5);
    %view(2)
    hold on
end
[Xr,Yr,Fr] = griddata(COORD(:,1),COORD(:,2),fieldp,Xp,Yp);
%vlevel=[0:.2:2];
%contourf(Xr,Yr,Fr,vlevel)
c2p=['r';'b';'g'];
if vel=='c'
    contour(Xr,Yr,Fr,20)
else
    contour(Xr,Yr,Fr,20,'LineColor',c2p(lay))
end

%colorbar
if isempty(data.hole_nodes1)
    hole=[];
else
    hole=convhull(COORD(data.hole_nodes1,:));
end
hold on
%vmM=max(vmM,max(max(Fd)));
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

