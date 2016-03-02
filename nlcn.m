%non linear constraints for all the layers
function [cons,ceq,dcons,dceq]=nlcn(dvs,G,data,nmin,nmax,pow,rmin,cdiv,rc,vfc)
%upload data 
nd=data.nd;
ne=data.N_ELEM;
nlay=data.nlay;
ceq=[];
dceq=[];
cons=[];
dcons=[];

%scale level set functiion
dv=dvs;
dv(1:(nlay*nd))=0.03*dvs(1:(nlay*nd));
dvf=G*dv;   %filter
dve=dvf((nd*nlay+1):end);
Area_elv=sparse(nlay,(ne+nd)*nlay);
At=sum(data.Area_el);
for lay=1:nlay
    Area_elv(lay,nlay*nd+(lay-1)*ne+(1:ne))=data.Area_el;
end
%volume fraction for each layer
vf=Area_elv*dvf/At-vfc;

%each layer
for lay=1:nlay
    dvl=dvf((lay-1)*nd+(1:nd));
    [consl,dconsl]=nlcnl(dvl,data,nmin,nmax,pow,rmin,lay,cdiv,rc);
    cons=[cons;consl];
    dcons=[dcons,dconsl];
end
%add volume fraction constraint
cons=[cons;vf];
dcons=[dcons,Area_elv'/At];
%filter & scale
dcons=G'*dcons;
dcons(1:(nd*nlay),:)=0.03*dcons(1:(nd*nlay),:); 

end
