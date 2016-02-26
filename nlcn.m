%non linear constraints for all the layers
function [cons,ceq,dcons,dceq]=nlcn(dvs,G,data,nmin,nmax,pow,rmin,cdiv,rc)
%scale level set functiion
dv=0.03*dvs;
nd=data.nd;
ceq=[];
dceq=[];
cons=[];
dcons=[];
dvf=G*dv;
%each layer
for lay=1:data.nlay
    dvl=dvf((lay-1)*nd+(1:nd));
    [consl,dconsl]=nlcnl(dvl,data,nmin,nmax,pow,rmin,lay,cdiv,rc);
    cons=[cons;consl];
    dcons=[dcons,dconsl];
end
dcons=0.03*G'*dcons;
end
