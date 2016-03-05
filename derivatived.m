%Sensitivity analysis function
function dtheta=derivatived(dv,WG,UG,data,ELEM_NODE,COORD,th,nmax,nlay,vel)
nd=data.nd;
ne=data.N_ELEM;
matC0=data.matC0;
dve=dv((nd*nlay+1):end);
dtheta=zeros(1,nlay*(nd+ne));
BKe=zeros(4,8);
NKe=zeros(2,8);
Ngauss=2;   %Number of Gauss points
[egv,wg] = GLTable(Ngauss);

%orthotropic material
c1111=matC0(1,1);
c2222=matC0(4,4);
c1122=matC0(1,4);
c2121=matC0(2,2);


%integration in the domain
for ele=1:ne
    %Index of element nodes vector field
    ienv=ELEM_NODE(ele,:);
    %index of nodes scalar field
    iens=ienv(2:2:end)/2;
    Xe=COORD(iens,:)';
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
            %kronecker
            BKe([1,3],1:2:end)=Be;
            BKe([2,4],2:2:end)=Be;
            %derivative of the displacement
            Du=BKe*UG(ienv,:);
            Dw=BKe*WG(ienv,:);
            
            %layer by layer
            for lay=1:nlay
                %angle by level set function
                %[angle, normnphi,dangle, ienbl, dnormnphi]=lsangle(x,data,ELEM_NODEb,dv((lay-1)*nd+(1:nd)),bs);
                dvl=dv((lay-1)*nd+(1:nd));
                %gradient level set function
                dphi=Be*dvl(iens);
                %normal of the gradient (diferentiable at 0)
                normnphi=sqrt((dphi'*dphi));
                %derivative of the gradient
                dnormnphi=dphi'*Be/normnphi;
                
                %angle with correct sign
                angle=atan2(-dphi(1),dphi(2));
                %derivative of the angle
                dangle=(-dphi(2)*Be(1,:)+dphi(1)*Be(2,:))/(dphi'*dphi);
                
                %index of element nodes
                ienl=(lay-1)*nd+iens;
                %rotated material matrix
                matC=[ c2121*(1-cos(4*angle))/2 + (c1111*(cos(2*angle)/2 + 1/2) - c1122*(cos(2*angle)/2 - 1/2))*(cos(2*angle)/2 + 1/2) - (c1122*(cos(2*angle)/2 + 1/2) - c2222*(cos(2*angle)/2 - 1/2))*(cos(2*angle)/2 - 1/2), (sin(2*angle)*(c1111 - c2222 + c1111*cos(2*angle) - 2*c1122*cos(2*angle) - 4*c2121*cos(2*angle) + c2222*cos(2*angle)))/4, (sin(2*angle)*(c1111 - c2222 + c1111*cos(2*angle) - 2*c1122*cos(2*angle) - 4*c2121*cos(2*angle) + c2222*cos(2*angle)))/4,                                                      c1111/8 + (3*c1122)/4 - c2121/2 + c2222/8 - (c1111*cos(4*angle))/8 + (c1122*cos(4*angle))/4 + (c2121*cos(4*angle))/2 - (c2222*cos(4*angle))/8;...
                    (sin(2*angle)*(c1111 - c2222 + c1111*cos(2*angle) - 2*c1122*cos(2*angle) - 4*c2121*cos(2*angle) + c2222*cos(2*angle)))/4,            c2121 + (c1111*sin(2*angle)^2)/4 - (c1122*sin(2*angle)^2)/2 - c2121*sin(2*angle)^2 + (c2222*sin(2*angle)^2)/4,            c2121 + (c1111*(1-cos(4*angle)))/8 - (c1122*(1-cos(4*angle)))/4 - c2121*(1-cos(4*angle))/2 + (c2222*(1-cos(4*angle)))/8,                                                                           (sin(2*angle)*(c1111 - c2222 - c1111*cos(2*angle) + 2*c1122*cos(2*angle) + 4*c2121*cos(2*angle) - c2222*cos(2*angle)))/4;...
                    (sin(2*angle)*(c1111 - c2222 + c1111*cos(2*angle) - 2*c1122*cos(2*angle) - 4*c2121*cos(2*angle) + c2222*cos(2*angle)))/4,            c2121 + (c1111*sin(2*angle)^2)/4 - (c1122*sin(2*angle)^2)/2 - c2121*sin(2*angle)^2 + (c2222*sin(2*angle)^2)/4,            c2121 + (c1111*(1-cos(4*angle)))/8 - (c1122*(1-cos(4*angle)))/4 - c2121*(1-cos(4*angle))/2 + (c2222*(1-cos(4*angle)))/8,                                                                           (sin(2*angle)*(c1111 - c2222 - c1111*cos(2*angle) + 2*c1122*cos(2*angle) + 4*c2121*cos(2*angle) - c2222*cos(2*angle)))/4;...
                    c1111/8 + (3*c1122)/4 - c2121/2 + c2222/8 - (c1111*cos(4*angle))/8 + (c1122*cos(4*angle))/4 + (c2121*cos(4*angle))/2 - (c2222*cos(4*angle))/8, (sin(2*angle)*(c1111 - c2222 - c1111*cos(2*angle) + 2*c1122*cos(2*angle) + 4*c2121*cos(2*angle) - c2222*cos(2*angle)))/4, (sin(2*angle)*(c1111 - c2222 - c1111*cos(2*angle) + 2*c1122*cos(2*angle) + 4*c2121*cos(2*angle) - c2222*cos(2*angle)))/4, c2121*(1-cos(4*angle))/2 + (c1111*(cos(2*angle)/2 - 1/2) - c1122*(cos(2*angle)/2 + 1/2))*(cos(2*angle)/2 - 1/2) - (c1122*(cos(2*angle)/2 - 1/2) - c2222*(cos(2*angle)/2 + 1/2))*(cos(2*angle)/2 + 1/2)];
                %derivative of the rotated material matrix respect to the angle
                dmatC=[                         -sin(2*angle)*(c1111 - c2222 + c1111*cos(2*angle) - 2*c1122*cos(2*angle) - 4*c2121*cos(2*angle) + c2222*cos(2*angle)), (c1111*cos(2*angle))/2 + (c1111*cos(4*angle))/2 - c1122*cos(4*angle) - 2*c2121*cos(4*angle) - (c2222*cos(2*angle))/2 + (c2222*cos(4*angle))/2, (c1111*cos(2*angle))/2 + (c1111*cos(4*angle))/2 - c1122*cos(4*angle) - 2*c2121*cos(4*angle) - (c2222*cos(2*angle))/2 + (c2222*cos(4*angle))/2,                                                                                          (sin(4*angle)*(c1111 - 2*c1122 - 4*c2121 + c2222))/2;...
                    (c1111*cos(2*angle))/2 + (c1111*cos(4*angle))/2 - c1122*cos(4*angle) - 2*c2121*cos(4*angle) - (c2222*cos(2*angle))/2 + (c2222*cos(4*angle))/2,                                                                                          (sin(4*angle)*(c1111 - 2*c1122 - 4*c2121 + c2222))/2,                                                                                          (sin(4*angle)*(c1111 - 2*c1122 - 4*c2121 + c2222))/2, (c1111*cos(2*angle))/2 - (c1111*cos(4*angle))/2 + c1122*cos(4*angle) + 2*c2121*cos(4*angle) - (c2222*cos(2*angle))/2 - (c2222*cos(4*angle))/2;...
                    (c1111*cos(2*angle))/2 + (c1111*cos(4*angle))/2 - c1122*cos(4*angle) - 2*c2121*cos(4*angle) - (c2222*cos(2*angle))/2 + (c2222*cos(4*angle))/2,                                                                                          (sin(4*angle)*(c1111 - 2*c1122 - 4*c2121 + c2222))/2,                                                                                          (sin(4*angle)*(c1111 - 2*c1122 - 4*c2121 + c2222))/2, (c1111*cos(2*angle))/2 - (c1111*cos(4*angle))/2 + c1122*cos(4*angle) + 2*c2121*cos(4*angle) - (c2222*cos(2*angle))/2 - (c2222*cos(4*angle))/2;...
                    (sin(4*angle)*(c1111 - 2*c1122 - 4*c2121 + c2222))/2, (c1111*cos(2*angle))/2 - (c1111*cos(4*angle))/2 + c1122*cos(4*angle) + 2*c2121*cos(4*angle) - (c2222*cos(2*angle))/2 - (c2222*cos(4*angle))/2, (c1111*cos(2*angle))/2 - (c1111*cos(4*angle))/2 + c1122*cos(4*angle) + 2*c2121*cos(4*angle) - (c2222*cos(2*angle))/2 - (c2222*cos(4*angle))/2,                          sin(2*angle)*(c1111 - c2222 - c1111*cos(2*angle) + 2*c1122*cos(2*angle) + 4*c2121*cos(2*angle) - c2222*cos(2*angle))];
                %if constant velocity, volume fractions are set to one
                if vel=='c'
                    %volume fraction
                    chi=normnphi/nmax;
                    %penalization
                    [rhol,drhol]=penal(chi);
                else
                    rhol=1;
                    drhol=0;
                end
                %volume fraction independent of level set
                chie=dve(ele+(lay-1)*ne);
                
                %material properties
                dmatC1=chie*rhol*dmatC/nlay;
                dmatC2=chie*drhol*matC/nmax/nlay;
                dmatC3=rhol*matC/nlay;
                
                %derivative
                dtheta(ienl)=dtheta(ienl)+wg(eit)*wg(nit)*Jdet*th*((Du'*dmatC2*Du/2-Dw'*dmatC2*Du)*dnormnphi+(Du'*dmatC1*Du/2-Dw'*dmatC1*Du)*dangle);
                dtheta(nd*nlay+ele+(lay-1)*ne)=dtheta(nd*nlay+ele+(lay-1)*ne)+wg(eit)*wg(nit)*Jdet*th*(Du'*dmatC3*Du/2-Dw'*dmatC3*Du);
                
            end
        end
    end
end
end
