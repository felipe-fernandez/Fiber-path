%Check sensitivity analysis with finite differences
function checkfindiff(obFUN,conFUN,data,numv)
%index of variable to check
iv=round((data.nd+data.N_ELEM)*data.nlay*rand(numv,1))+1;

%Initial guess
%dv=0.9*(data.Yc+15/1000);%rand(size(data.Yb));%data.Yb,yH-data.Yb);%min(min(data.Yb,data.Xb),min(yH-data.Yb,xL-data.Xb));
dv=sqrt(data.Yc.^2+data.Xc.^2);
dv=repmat(dv,data.nlay,1);
dv=[dv;0.5*ones(data.N_ELEM*data.nlay,1)];
%initial fea to normalize
[theta,dtheta]=obFUN(dv);
[cons,ceq,dcons,dceq]=conFUN(dv);

%finite differences
ddiff=1e-4;
for ivar=1:numv
    alphad=dv;
    var=iv(ivar);
    alphad(var)=alphad(var)+ddiff;
    %finite element and post processing function
    [thetad,dthetad]=obFUN(alphad);
    [consd,ceqd,dconsd,dceqd]=conFUN(alphad);
    %derivative
    dtheta_fd(ivar)=(thetad-theta)/ddiff;
    dcons_fd(:,ivar)=(consd-cons)/ddiff;
end

display('--------------------------------------------------')
display('Check sensitivity analysis with finite differences')
display('--------------------------------------------------')
display(num2str([dtheta_fd' dtheta(iv)'],'%14.5e'))
display(['error norm =' num2str(norm(abs(dtheta_fd'-dtheta(iv)')),'%14.5e')])
for it=1:size(dcons,2)
    display('--------------------------------------------------')
    display(num2str([dcons_fd(it,:)' dcons(iv,it)],'%14.5e'))
    display(['error norm =' num2str(norm(abs(dcons_fd(it,:)'-dcons(iv,it))),'%14.5e')])
end
display('--------------------------------------------------')
end