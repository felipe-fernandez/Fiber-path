%Finite element function
function [theta,dtheta]=feafun(dvs,G,data,UG0,FG,th,c0,nmax)

ELEM_NODE=data.ELEM_NODE;
COORD=data.COORD;
nameplot=data.nameplot;
nlay=data.nlay;

dv=0.03*dvs;
nd=data.nd; %number of design variables

%FINITE ELEMENT ANALYSIS
[UG,KG]=FEAsolver(data,ELEM_NODE,UG0,FG,G*dv,COORD,th,nmax,nlay);

%Constitutive matrix
matC0=data.matC0;
stressana(dv,data,ELEM_NODE,matC0,COORD,UG,nlay,nmax);
%adjoint problem
[WG,theta]=adjointFEA(dv,UG,data,ELEM_NODE,COORD,matC0,th,KG,nlay,nmax);
%self adjoint problem
%WG=UG;
%theta=FG'*UG/2;

%derivative
dtheta=derivatived(dv,WG,UG,data,ELEM_NODE,COORD,matC0,th,nmax,nlay);
%normalize objective
theta=theta/c0;
%normalize and filter chain rule
dtheta=0.03*dtheta*G/c0;
%filter chain rule
dvf=G*dv;

%plot each layer
for lay=1:nlay
    plotfun(COORD,dvf((lay-1)*nd+(1:nd)),data,nameplot,lay);
    hold off
end

for lay=1:nlay
    plotfun(COORD,dvf((lay-1)*nd+(1:nd)),data,nameplot,nlay+1);
end
hold off
end
