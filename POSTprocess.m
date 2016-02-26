%objective and adjoint problem for topological derivative
function [theta,dtheta]=POSTprocess(UG,KG,COORD,ELEM_NODE,data,th,dv,FG,nmax,nlay)
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
end
