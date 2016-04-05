function [ptof,change]=replaceboundary(ptoi,ptof,l,ptob1,ptob2,mb)
change=0;
%slope
mp=(ptof(2)-ptoi(2))/(ptof(1)-ptoi(1));

%box points
if abs(mb)>0
    unitp=[1,-1/mb]/(norm([1,-1/mb]));
else
    unitp=[0,1];
end
boxp=[ptob1-l*unitp;ptob2-l*unitp;ptob2+l*unitp;ptob1+l*unitp];

%check if points are inside the box
in = inpolygon([ptoi(1);ptof(1)],[ptoi(2);ptof(2)],boxp(:,1),boxp(:,2));

%if it is inside then it is close to our boundary
if in(1)==1 || in(2)==1
    
    %mp differ mb
    x=(ptoi(2)-ptob1(2)-ptoi(1)*mp+ptob1(1)*mb)/(mb-mp);
    y=ptoi(2)+(x-ptoi(1))*mp;
    
    %replace
    ptof=[x;y];
    change=1;
    %xpto=linspace(ptoi(1),ptof(1),100);
    %ypto=ptoi(2)+(xpto-ptoi(1))*mp;
    %xptob=linspace(ptob1(1),ptob2(1),100);
    %yptob=ptob1(2)+(xptob-ptob1(1))*mb;
    
    %     figure()
    %     plot(x,y,'ro',xpto,ypto,'b-',xptob,yptob,'g-',boxp(:,1),boxp(:,2),'*');
end
end