%non linear constraints for all the layers
function [cons,ceq,dcons,dceq,vfrac]=nlcn(dvs,G,data,nmin,nmax,pow,rmin,cdiv,rc)
ELEM_NODE=data.ELEM_NODE;
COORD=data.COORD;
nlay=data.nlay;

dv=0.03*dvs;
nd=data.nd;
ceq=[];
dceq=[];
cons=[];
dcons=[];
dvf=G*dv;
vfrac=0;
%each layer
for lay=1:nlay
    dvl=dvf((lay-1)*nd+(1:nd));
    [consl,dconsl,vfracl]=nlcnl(dvl,data,ELEM_NODE,nmin,nmax,pow,rmin,lay,nlay,COORD,cdiv,rc);
    cons=[cons;consl];
    dcons=[dcons,dconsl];
    vfrac=vfrac+vfracl/nlay;
end
dcons=0.03*G'*dcons;
end
