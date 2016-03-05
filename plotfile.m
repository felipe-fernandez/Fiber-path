%Function that reads file and plots 
function plotfile(nlay, rc, efilter, vel)

nameplot=['Fig' vel num2str(nlay) rc num2str(round(efilter*1000))]; %Fig nlay constraints filter
load([nameplot '.mat'])

% nd=data.nd;
% ne=data.N_ELEM;
% nlay=data.nlay;
% COORD=data.COORD;
% 
% dvs=dvo;
% dv=dvs;
% dv(1:(nlay*nd))=0.03*dvs(1:(nlay*nd)); %scale design variables
% dvf=G*dv;       %filter design variables
% 
%     dve=dvf((nd*nlay+1):end);
%     dva=zeros(ne,1);
%     %plot each layer
%     for lay=1:nlay
%         plotfun(COORD,dvf((lay-1)*nd+(1:nd)),dve((lay-1)*ne+(1:ne)),data,nameplot,lay,lay);
%         hold off
%         dva=dva+dve((lay-1)*ne+(1:ne))/nlay;
%     end
%     
%     for lay=1:nlay
%         plotfun(COORD,dvf((lay-1)*nd+(1:nd)),dva,data,nameplot,nlay+1,lay);
%     end
%     hold off
    
[theta,dtheta]=feafun(dvo,G,data,UG0,FG,th,c0,nmax,1);
[cons,ceq,dcons,dceq]=nlcn(dvo,G,data,nmin,nmax,pow,rmin,cdiv,rc,vfc);
disp([nameplot '  g=' num2str(nmin) '/' num2str(nmax) '  filter=' num2str(efilter) '  rmin=' num2str(rmin)  '  cdiv='  num2str(cdiv) '  c=' num2str(theta,'%1.8f') ])
end
