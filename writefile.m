function writefile(C,num,nameplot,thic)
%read contour matrix and obtain points
lastce=size(C,2);    % last entry of the Cmatrix
cen=1;
numt=0;
D=[];

%box thickness
l=5e-3;
%top line horizontal line from (x1,y1) to (x2,y2)
ptob1=[-15e-3,15e-3-thic/2];
ptob2=[15e-3,15e-3-thic/2];
%slope
mb12=(ptob2(2)-ptob1(2))/(ptob2(1)-ptob1(1));
%bottom line horizontal line from (x1,y1) to (x2,y2)
ptob3=[-15e-3,-15e-3+thic/2];
ptob4=[15e-3,-15e-3+thic/2];
%slope
mb34=(ptob4(2)-ptob3(2))/(ptob4(1)-ptob3(1));
%center of left circle
xo1=[-15e-3;0];
xo2=[15e-3;0];
Rbig=(15e-3)-thic/2;
rsma=(5e-3)+thic/2;

plopo=1;
%loop over toolpahts
while cen<lastce
    %number of points
    numxy=C(2,cen);
    if numxy<10
        %no line if it is made of just three points
    else
        %line with more than one point
        XYP=C(:,cen+(1:numxy));
        
        %first points
        ptoi=XYP(:,2);
        ptof=XYP(:,1);
        %replace boundary top
        [XYP(:,1),change]=replaceboundary(ptoi,ptof,l,ptob1,ptob2,mb12);
        if change==0
            %bottom
            [XYP(:,1),change]=replaceboundary(ptoi,ptof,l,ptob3,ptob4,mb34);
            if change==0
                %left circle
                [XYP(:,1),change]=replaceoutcircle(ptoi,ptof,xo1,Rbig,l,'l');
                if change==0
                    %right circle
                    [XYP(:,1),change]=replaceoutcircle(ptoi,ptof,xo2,Rbig,l,'r');
                    if change==0
                        %right circle
                        [XYP(:,1),change]=replaceoutcircle(ptoi,ptof,xo1,rsma,l,'o');
                        if change==0
                            %right circle
                            [XYP(:,1),change]=replaceoutcircle(ptoi,ptof,xo2,rsma,l,'o');
                        end
                    end
                end
            end
        end
        
                
        %last points
        ptoi=XYP(:,end-1);
        ptof=XYP(:,end);
        %replace boundary
        [XYP(:,end),change]=replaceboundary(ptoi,ptof,l,ptob1,ptob2,mb12);
        if change==0
            [XYP(:,end),change]=replaceboundary(ptoi,ptof,l,ptob3,ptob4,mb34);
            if change==0
                [XYP(:,end),change]=replaceoutcircle(ptoi,ptof,xo1,Rbig,l,'l');
                if change==0
                    [XYP(:,end),change]=replaceoutcircle(ptoi,ptof,xo2,Rbig,l,'r');
                    if change==0
                        [XYP(:,end),change]=replaceoutcircle(ptoi,ptof,xo1,rsma,l,'o');
                        if change==0
                            [XYP(:,end),change]=replaceoutcircle(ptoi,ptof,xo2,rsma,l,'o');
                        end
                    end
                end
            end
        end
        
        
        %distance between points       
        dx=XYP(1,2:end)-XYP(1,1:(end-1));
        dy=XYP(2,2:end)-XYP(2,1:(end-1));
        dz=sqrt(dx.^2+dy.^2);
        
        %check distance between points
        auxvec=2:(numxy-1);
        
        %aproved points
        apoi=[1 auxvec(dz(1:(end-1))>5e-4 & [ones(1,numxy-3) dz(end)>5e-4]) numxy];
        
        %newpoints that are at least 5e-4 apart
        nXYP=[XYP(:,apoi)];
        D=[D [0;size(nXYP,2)] nXYP];
        %number of toolpaths counter
        numt=numt+1;
    end
    %next line
    cen=numxy+cen+1;
end
if plopo==1
    figure(num+3)
end
cen=1;
%write in file with the same nameplot
namefilep=[nameplot num2str(num)];
fio=fopen([namefilep '.txt'],'w');
fprintf(fio,'%s\n',['# test file for ' namefilep]);
%total number of toolpaths
fprintf(fio,'%d\n',numt);
lastde=size(D,2);    % last entry of the Cmatrix
%loop over toolpahts
while cen<lastde
    %number of points
    numxy=D(2,cen);
    %line with more than one point
    XYP=D(:,cen+(1:numxy));
    fprintf(fio,'%d\n',numxy);
    fprintf(fio,'%+9.6f %+9.6f\n',XYP);
    if plopo==1
        plot(XYP(1,:),XYP(2,:),'-+');%,'Markersize',3,'Markerfacecolor','b')
        
        hold on
    end
    %next line
    cen=numxy+cen+1;
end
hold off
%close file
fclose(fio);

if plopo==1
    axis equal
    axis tight
    set(gcf,'paperunits','centimeters')
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf,'papersize',[16,16])
    set(gcf,'paperposition',[0,0,16,16])
    saveas(gcf,[nameplot num2str(num) '.eps']);
end
end