%Finite element function
function [theta,dtheta]=feafun(dvs,G,data,UG0,FG,th,c0,nmax,plv,vel)

%upload data
ELEM_NODE=data.ELEM_NODE;
COORD=data.COORD;
nameplot=data.nameplot;
nlay=data.nlay;
nd=data.nd; %number of design variables
ne=data.N_ELEM;

%if constant velocity, volume fractions are set to one
if vel=='c'
    dvs((nlay*nd+1):(nd+ne)*nlay)=1;
end
dv=dvs;
dv(1:(nlay*nd))=0.03*dvs(1:(nlay*nd)); %scale design variables
dvf=G*dv;       %filter design variables

%FINITE ELEMENT ANALYSIS
[UG,KG]=FEAsolver(data,ELEM_NODE,UG0,FG,dvf,COORD,th,nmax,nlay,vel);

%stressana(dvf,data,ELEM_NODE,matC0,COORD,UG,nlay,nmax);
%adjoint problem
[WG,theta]=adjointFEA(dvf,UG,data,ELEM_NODE,COORD,th,KG,nlay,nmax,vel);
%self adjoint problem
%WG=UG;
%theta=FG'*UG/2;

%derivative
dtheta=derivatived(dvf,WG,UG,data,ELEM_NODE,COORD,th,nmax,nlay,vel);
%normalize objective
theta=theta/c0;
%normalize and filter chain rule
dtheta=dtheta*G/c0;
dtheta(1:(nd*nlay))=0.03*dtheta(1:(nd*nlay));

%if constant velocity the design variables are just the level set function
if vel=='c'
    dtheta=dtheta(1:(nd*nlay));
end

if plv==1
    dve=dvf((nd*nlay+1):end);
    dva=zeros(ne,1);
    %plot each layer
    for lay=1:nlay
        plotfun(COORD,dvf((lay-1)*nd+(1:nd)),dve((lay-1)*ne+(1:ne)),data,nameplot,lay,lay,vel);
        hold off
        dva=dva+dve((lay-1)*ne+(1:ne))/nlay;
    end
    
    for lay=1:nlay
        plotfun(COORD,dvf((lay-1)*nd+(1:nd)),dva,data,nameplot,nlay+1,lay,vel);
    end
    hold off
end
end
