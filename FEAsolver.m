%Finite Element Analysis function solve the solid 2D problem
function [UG,KG]=FEAsolver(data,ELEM_NODE,UG,FG,dv,COORD,th,nmax,nlay)
%Loop over elements to assembly the stiffness matrix
KG=sparse(data.N_NODE,data.N_NODE); %Global Stiffness matrix
Ngauss=2;   %Number of Gauss points
[egv,wg] = GLTable(Ngauss);
ngv=egv;
matC0=data.matC0;

%orthotropic material
c1111=matC0(1,1);
c2222=matC0(4,4);
c1122=matC0(1,4);
c2121=matC0(2,2);

BKe=zeros(4,8);
nd=data.nd;

%density mesh
for ele=1:data.N_ELEM
    %Index of element nodes vector field
    ienv=ELEM_NODE(ele,:);
    %index of nodes scalar field
    ine=ienv(2:2:end)/2;
    %coordinates of the points
    Xe=COORD(ine,:)';
    %Gauss quadrature integration
    for eit=1:Ngauss
        eg=egv(eit);
        for nit=1:Ngauss
            ng=ngv(nit);
            
            %Derivative of the shape function
            DNe=1/4*[-(1-ng) 1-ng 1+ng -(1+ng);...
                -(1-eg) -(1+eg) 1+eg 1-eg]';
            %Jacobian of the element e at location e,n
            J=(Xe*DNe);
            Jdet=det(J);
            %Derivative of the shape function
            Be=(DNe/J)';
            matC=zeros(4,4);
            %layer by layer
            for lay=1:nlay
                %angle by level set function 1
                %[angle,normnphi]=lsangle(x,data,ELEM_NODEb,dv((lay-1)*nd+(1:nd)),bs);
                dvl=dv((lay-1)*nd+(1:nd));
                %gradient level set function
                dphi=Be*dvl(ine);
                %normal of the gradient (diferentiable at 0)
                normnphi=sqrt((dphi'*dphi));
                
                if normnphi<1e-3
                    'Error: Flat level set surface'
                end
                
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
            
            %Assembly
            %Stifness for element in domain
            KG(ienv,ienv)=KG(ienv,ienv)+wg(eit)*wg(nit)*BKe'*matC*BKe*Jdet*th;
        end
    end
end
%Solve the displacement in the freenodes FEA(freedofs,:)
UG(data.free,:)=KG(data.free,data.free)\(FG(data.free,:)...
    -KG(data.free,data.pres)*UG(data.pres,:));
FG(data.pres)=KG(data.pres,data.free)*UG(data.free,:)+...
    KG(data.pres,data.pres)*UG(data.pres,:);

end
