function [ptof,change]=replaceoutcircle(ptoi,ptof,xo,r,l,pos)
 change=0;
%delta x
dx=ptof-xo;
%radius of the point respect to xo
rp=norm(dx);

%inside region
if abs(rp-r)<l && ((pos=='l' && ptof(1)<xo(1)) || (pos=='r' && ptof(1)>xo(1)) || pos=='o')
    
    %slope
    mp=(ptof(2)-ptoi(2))/(ptof(1)-ptoi(1));
    b=(ptoi(2)+ptof(2)-mp*(ptoi(1)+ptof(1)))/2;
    %cross
    xcf1=(xo(1) - b*mp + mp*xo(2) - (- b^2 - 2*b*mp*xo(1) + 2*b*xo(2) + mp^2*r^2 - mp^2*xo(1)^2 + 2*mp*xo(1)*xo(2) + r^2 - xo(2)^2)^(1/2))/(mp^2 + 1);
    xcf2=(xo(1) - b*mp + mp*xo(2) + (- b^2 - 2*b*mp*xo(1) + 2*b*xo(2) + mp^2*r^2 - mp^2*xo(1)^2 + 2*mp*xo(1)*xo(2) + r^2 - xo(2)^2)^(1/2))/(mp^2 + 1);
    
    if isreal(xcf1) && isreal(xcf2)
    
    pcf1=[xcf1;mp*xcf1+b];
    pcf2=[xcf2;mp*xcf2+b];
     %only one solution
    if norm(pcf1-ptof)> norm(pcf2-ptof)
        ptof=pcf2;
    else
        ptof=pcf1;
    end
    change=1;
    end
end

end