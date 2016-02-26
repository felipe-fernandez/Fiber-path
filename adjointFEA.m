
%adjoint problem
function [WG,theta]=adjointFEA(dv,UG,data,ELEM_NODE,COORD,matC0,th,KG,nlay,nmax)
Ngauss=2;   %Number of Gauss points
[egv,wg] = GLTable(Ngauss);
theta=0;
nd=data.nd;
FGw=zeros(data.N_NODE,1); %Global Stiffness matrix

BKe=zeros(4,8);
NKe=zeros(2,8);

%orthotropic material
c1111=matC0(1,1);
c2222=matC0(4,4);
c1122=matC0(1,4);
c2121=matC0(2,2);

%integration in the domain
for ele=1:data.N_ELEM
    %Index of element nodes vector field
    ienv=ELEM_NODE(ele,:);
    %index of nodes scalar field
    ine=ienv(2:2:end)/2;
    Xe=COORD(ine,:)';
    %Gauss quadrature integration
    for eit=1:Ngauss
        eg=egv(eit);
        for nit=1:Ngauss
            ng=egv(nit);
            %Shape functions
            Ne=1/4*(1+eg*[-1 1 1 -1]).*(1+ng*[-1 -1 1 1]);
            %Derivative of the shape function
            DNe=1/4*[-(1-ng) 1-ng 1+ng -(1+ng);...
                -(1-eg) -(1+eg) 1+eg 1-eg]';
            Jt=Xe*DNe;
            Jdet=det(Jt);
            %displacement
            NKe(1,1:2:end)=Ne;
            NKe(2,2:2:end)=Ne;
            %Derivative of the shape function
            Be=(DNe/Jt)';
            matC=zeros(4,4);
            %layer by layer
            for lay=1:nlay
                %nodal level set function for layer
                dvl=dv((lay-1)*nd+(1:nd));
                %gradient level set function
                dphi=Be*dvl(ine);
                %normal of the gradient (diferentiable at 0)
                normnphi=sqrt((dphi'*dphi));
                
                %angle with correct sign
                angle=atan2(-dphi(1),dphi(2));
                
                %material matrix
                matCa=[ c2121*(1-cos(4*angle))/2 + (c1111*(cos(2*angle)/2 + 1/2) - c1122*(cos(2*angle)/2 - 1/2))*(cos(2*angle)/2 + 1/2) - (c1122*(cos(2*angle)/2 + 1/2) - c2222*(cos(2*angle)/2 - 1/2))*(cos(2*angle)/2 - 1/2), (sin(2*angle)*(c1111 - c2222 + c1111*cos(2*angle) - 2*c1122*cos(2*angle) - 4*c2121*cos(2*angle) + c2222*cos(2*angle)))/4, (sin(2*angle)*(c1111 - c2222 + c1111*cos(2*angle) - 2*c1122*cos(2*angle) - 4*c2121*cos(2*angle) + c2222*cos(2*angle)))/4,                                                      c1111/8 + (3*c1122)/4 - c2121/2 + c2222/8 - (c1111*cos(4*angle))/8 + (c1122*cos(4*angle))/4 + (c2121*cos(4*angle))/2 - (c2222*cos(4*angle))/8;...
                    (sin(2*angle)*(c1111 - c2222 + c1111*cos(2*angle) - 2*c1122*cos(2*angle) - 4*c2121*cos(2*angle) + c2222*cos(2*angle)))/4,            c2121 + (c1111*sin(2*angle)^2)/4 - (c1122*sin(2*angle)^2)/2 - c2121*sin(2*angle)^2 + (c2222*sin(2*angle)^2)/4,            c2121 + (c1111*(1-cos(4*angle)))/8 - (c1122*(1-cos(4*angle)))/4 - c2121*(1-cos(4*angle))/2 + (c2222*(1-cos(4*angle)))/8,                                                                           (sin(2*angle)*(c1111 - c2222 - c1111*cos(2*angle) + 2*c1122*cos(2*angle) + 4*c2121*cos(2*angle) - c2222*cos(2*angle)))/4;...
                    (sin(2*angle)*(c1111 - c2222 + c1111*cos(2*angle) - 2*c1122*cos(2*angle) - 4*c2121*cos(2*angle) + c2222*cos(2*angle)))/4,            c2121 + (c1111*sin(2*angle)^2)/4 - (c1122*sin(2*angle)^2)/2 - c2121*sin(2*angle)^2 + (c2222*sin(2*angle)^2)/4,            c2121 + (c1111*(1-cos(4*angle)))/8 - (c1122*(1-cos(4*angle)))/4 - c2121*(1-cos(4*angle))/2 + (c2222*(1-cos(4*angle)))/8,                                                                           (sin(2*angle)*(c1111 - c2222 - c1111*cos(2*angle) + 2*c1122*cos(2*angle) + 4*c2121*cos(2*angle) - c2222*cos(2*angle)))/4;...
                    c1111/8 + (3*c1122)/4 - c2121/2 + c2222/8 - (c1111*cos(4*angle))/8 + (c1122*cos(4*angle))/4 + (c2121*cos(4*angle))/2 - (c2222*cos(4*angle))/8, (sin(2*angle)*(c1111 - c2222 - c1111*cos(2*angle) + 2*c1122*cos(2*angle) + 4*c2121*cos(2*angle) - c2222*cos(2*angle)))/4, (sin(2*angle)*(c1111 - c2222 - c1111*cos(2*angle) + 2*c1122*cos(2*angle) + 4*c2121*cos(2*angle) - c2222*cos(2*angle)))/4, c2121*(1-cos(4*angle))/2 + (c1111*(cos(2*angle)/2 - 1/2) - c1122*(cos(2*angle)/2 + 1/2))*(cos(2*angle)/2 - 1/2) - (c1122*(cos(2*angle)/2 - 1/2) - c2222*(cos(2*angle)/2 + 1/2))*(cos(2*angle)/2 + 1/2)];
                %volume fraction
                chi=normnphi/nmax;
                %penalization
                rho=penal(chi);
                %material properties
                matC=matC+rho*matCa/nlay;
            end
            
            %kronecker
            BKe([1,3],1:2:end)=Be;
            BKe([2,4],2:2:end)=Be;
            %stress
            Du=BKe*UG(ienv,:);
            S=matC*Du;
            
            
            %compliance
            psi=S'*Du/2;
            npsi=BKe'*(matC'+matC)*Du/2;
            
            %assembly
            FGw(ienv,:)=FGw(ienv,:)+wg(eit)*wg(nit)*npsi*Jdet*th;
            %function
            theta=theta+wg(eit)*wg(nit)*psi*Jdet*th;
        end
    end
end
%Solve the displacement in the freenodes Adjoint problem
WG=zeros(data.N_NODE,1);
WG(data.free,:)=KG(data.free,data.free)\(FGw(data.free,:));
%FGw(data.pres,:)=KG(data.pres,data.free)*WG(data.free,:);
end
