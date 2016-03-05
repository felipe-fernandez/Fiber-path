% initialization function
% calls:    read file function to upload mesh
%           constitutive model function
%           cone filter function
% input: 
%   filename    - txt file with mesh
%   nlay        - number of layers
%   loadv       - nodal force
%   th          - thickness of the layer
%   rc          - requested constraints
%   efilter     - cone filter radius
%   fpres       - pressure force module
% output: described below
function [data,UG0,FG,G]=initializationf(filename,nlay,loadv,th,rc,efilter,fpres,vel)
%read file
data=readfile(filename);
%Constitutive matrix
data.matC0=consti();
%name of the plot and archives
data.nameplot=['Fig' vel num2str(nlay) rc num2str(round(efilter*1000))]; %Fig nlay constraints filter
data.At=pi*.015^2+.03^2-2*pi*.005^2;
data.nlay=nlay;
%cone filter
[G,data]=cfilter(data,efilter);

% Structure of the problem
FG=sparse(data.N_NODE,1);     %Global Applied Force Vector
UG0=sparse(data.N_NODE,1);     %Global Displacement Vector

%load at the node
FG(2*data.f_nodes1-1)=-loadv;
FG(2*data.f_nodes2-1)=loadv;

%Loop over elements to assembly the stiffness matrix
Ngauss=2;   %Number of Gauss points
[egv,wg] = GLTable(Ngauss);
ngv=egv;
NKe=zeros(2,8);

hside=0;
%elements and nodes to apply force
f_nodes=union(data.f_nodes1,data.f_nodes2);
f_ele=union(data.f_ele1,data.f_ele2)';
xplot=[];
forcp=[];

for ele=f_ele
    %Index of element nodes vector field
    ienv=data.ELEM_NODE(ele,:);
    %index of nodes scalar field
    ine=ienv(2:2:end)/2;
    %side recognition
    side=ismember(ine,f_nodes,'legacy');
    
    %coordinates of the points
    Xe=data.COORD(ine,:)';
    
    if side(2)==1 && side(3)==1
        %Gauss quadrature integration
        eg=1;
        for nit=1:Ngauss
            ng=ngv(nit);
            %Shape functions
            Ne=1/4*(1+eg*[-1 1 1 -1]).*(1+ng*[-1 -1 1 1]);
            %Derivative of the shape function
            DNe=1/4*[-(1-ng) 1-ng 1+ng -(1+ng);...
                -(1-eg) -(1+eg) 1+eg 1-eg]';
            %point
            xp=Ne*Xe';
            xplot=[xplot;xp];
            forcp=[forcp;xp(2) ((2.5/1000)^2-xp(2)^2)/((2.5/1000)^2)];
            forcem=((2.5/1000)^2-xp(2)^2)*sign(xp(1))*fpres/((2.5/1000)^2);
            %Jacobian of the element e at location e,n
            J=(Xe*DNe)';
            %dy
            dJ=sqrt(J(2,:)*J(2,:)');
            %kronecker
            NKe(1,1:2:end)=Ne;
            NKe(2,2:2:end)=Ne;
            %Assembly
            
            %force vector
            FG(ienv)=FG(ienv)+wg(nit)*NKe'*forcem*dJ*th;
            hside=hside+wg(nit)*dJ;
         end
    elseif side(2)==1 && side(1)==1
        %Gauss quadrature integration
        ng=-1;
        for nit=1:Ngauss
            eg=ngv(nit);
            %Shape functions
            Ne=1/4*(1+eg*[-1 1 1 -1]).*(1+ng*[-1 -1 1 1]);
            %point
            xp=Ne*Xe';
            xplot=[xplot;xp];
            forcp=[forcp;xp(2) ((2.5/1000)^2-xp(2)^2)/((2.5/1000)^2)];
            forcem=((2.5/1000)^2-xp(2)^2)*sign(xp(1))*fpres/((2.5/1000)^2);
            
            %Derivative of the shape function
            DNe=1/4*[-(1-ng) 1-ng 1+ng -(1+ng);...
                -(1-eg) -(1+eg) 1+eg 1-eg]';
            %Jacobian of the element e at location e,n
            J=(Xe*DNe)';
            %dy
            dJ=sqrt(J(1,:)*J(1,:)');
            %kronecker
            NKe(1,1:2:end)=Ne;
            NKe(2,2:2:end)=Ne;
            %Assembly
            %force vector
            FG(ienv)=FG(ienv)+wg(nit)*NKe'*forcem*dJ*th;
            hside=hside+wg(nit)*dJ;
        end
    elseif side(1)==1 && side(4)==1
        %Gauss quadrature integration
        eg=-1;
        for nit=1:Ngauss
            ng=ngv(nit);
            %Shape functions
            Ne=1/4*(1+eg*[-1 1 1 -1]).*(1+ng*[-1 -1 1 1]);
            %point
            xp=Ne*Xe';
            xplot=[xplot;xp];
            forcp=[forcp;xp(2) ((2.5/1000)^2-xp(2)^2)/((2.5/1000)^2)];
            forcem=((2.5/1000)^2-xp(2)^2)*sign(xp(1))*fpres/((2.5/1000)^2);
            
            %Derivative of the shape function
            DNe=1/4*[-(1-ng) 1-ng 1+ng -(1+ng);...
                -(1-eg) -(1+eg) 1+eg 1-eg]';
            %Jacobian of the element e at location e,n
            J=(Xe*DNe)';
            %dy
            dJ=sqrt(J(2,:)*J(2,:)');
            %kronecker
            NKe(1,1:2:end)=Ne;
            NKe(2,2:2:end)=Ne;
            %Assembly
            %force vector
            FG(ienv)=FG(ienv)+wg(nit)*NKe'*forcem*dJ*th;
            hside=hside+wg(nit)*dJ;
        end
    elseif side(3)==1 && side(4)==1
        %Gauss quadrature integration
        ng=1;
        for nit=1:Ngauss
            eg=ngv(nit);
            %Shape functions
            Ne=1/4*(1+eg*[-1 1 1 -1]).*(1+ng*[-1 -1 1 1]);
            %point
            xp=Ne*Xe';
            xplot=[xplot;xp];
            forcp=[forcp;xp(2) ((2.5/1000)^2-xp(2)^2)/((2.5/1000)^2)];
            forcem=((2.5/1000)^2-xp(2)^2)*sign(xp(1))*fpres/((2.5/1000)^2);
            
            %Derivative of the shape function
            DNe=1/4*[-(1-ng) 1-ng 1+ng -(1+ng);...
                -(1-eg) -(1+eg) 1+eg 1-eg]';
            %Jacobian of the element e at location e,n
            J=(Xe*DNe)';
            %dy
            dJ=sqrt(J(1,:)*J(1,:)');
            %kronecker
            NKe(1,1:2:end)=Ne;
            NKe(2,2:2:end)=Ne;
            %Assembly
            %force vector
            FG(ienv)=FG(ienv)+wg(nit)*NKe'*forcem*dJ*th;
            hside=hside+wg(nit)*dJ;
        end
    else
        side
    end
end
%plot(forcp(:,1),forcp(:,2))
end