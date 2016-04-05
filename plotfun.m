%plot scalar nodal values fieldp
function plotfun(COORD,fieldp,dene,data,nameplot,num,lay,vel,nmax,prob)
figure(num)
set(num,'Position',[20+400*(num-1) 20 400 400]);
vmM=max(fieldp);
%thickness of the bead
thic=6e-4;

xplot=linspace(min(data.Xc)+thic/2,max(data.Xc)-thic/2,120);
yplot=linspace(min(data.Yc)+thic/2,max(data.Yc)-thic/2,120);
[Xp,Yp]=meshgrid(xplot,yplot);
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
%contourf(Xr,Yr,Fr,vlevel)
limi=max(max(abs(Fr)));

%no holes
if prob=='c'
    radsmal=(5e-3)+thic/2;
    Radbig=(15e-3)-thic/2;
    Fr(sqrt((Xr+15e-3).^2+Yr.^2)<radsmal)=NaN;
    Fr(sqrt((Xr-15e-3).^2+Yr.^2)<radsmal)=NaN;
    Fr(sqrt((Xr+15e-3).^2+Yr.^2)>Radbig & Xr<-15e-3)=NaN;
    Fr(sqrt((Xr-15e-3).^2+Yr.^2)>Radbig & Xr>15e-3)=NaN;
end

clev=(-limi-nmax*thic):nmax*thic:(limi+nmax*thic);
c2p=['r';'b';'g'];
if vel=='c'
    [C,h]=contour(Xr,Yr,Fr,clev,'LineWidth',1.25);%,'LineColor',c2p(lay))
else
    [C,h]=contour(Xr,Yr,Fr,clev,'LineColor',c2p(lay));
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

% %write file
% if num~=(data.nlay+1)
%     writefile(C,num,nameplot,thic)
% end
end

