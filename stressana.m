%Nodal stress components averaging for nodes and plot of the vonMises
%stress
function stressana(dv,data,ELEM_NODE,matC0,COORD,UG,nlay,nmax)
nd=data.nd;
ne=data.N_ELEM;
dve=dv((nd*nlay+1):end);
BKe=zeros(4,8);
%evaluate at the nodes
Ngauss=1;   %Number of Gauss points
[ep,wg] = GLTable(Ngauss);
np=ep;
NODE_STRESS=zeros(4,nd);
NODES_NEL=zeros(4,nd);
%orthotropic material
c1111=matC0(1,1);
c2222=matC0(4,4);
c1122=matC0(1,4);
c2121=matC0(2,2);

for ele=1:ne
    %Index ol nodes for element el
    inev=ELEM_NODE(ele,:);
    %element matrix giving the coordinates of the element
    Xe=COORD(round(inev(1:2:end)/2),:)';
    %Stifness for element in domain 1
    ine=inev(2:2:end)/2;
            
    %Gauss quadrature integration
    for ptnodee=1:Ngauss
        eg=ep(ptnodee);
        for ptnoden=1:Ngauss
            ng=np(ptnoden);
            %Shape functions
            Ne=1/4*(1+eg*[-1 1 1 -1]).*(1+ng*[-1 -1 1 1]);
            %Derivative of the shape function
            DNe=1/4*[-(1-ng) 1-ng 1+ng -(1+ng);...
                -(1-eg) -(1+eg) 1+eg 1-eg]';
            %Jacobian of the element e at location e,n
            J=Xe*DNe;
            Jdet=det(J);
            %Derivative of the shape function
            Be=DNe/J;
            %kronecker
            BKe([1,3],1:2:end)=Be';
            BKe([2,4],2:2:end)=Be';
                        
            matC=zeros(4,4);
            %layer by layer
            for lay=1:nlay
                %nodal level set function for layer
                dvl=dv((lay-1)*nd+(1:nd));
                %gradient level set function
                dphi=Be'*dvl(ine);
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
                chil=normnphi/nmax;
                %penalization
                rhol=penal(chil);
                %volume fraction independent of level set
                chie=dve(ele+(lay-1)*ne);
                %material properties
                matC=matC+chie*rhol*matCa/nlay;
            end

            %stress
            Sv=matC*BKe*UG(inev);
            %nodal averaging
            NODE_STRESS(:,ine)=...
                NODE_STRESS(:,ine)+[Sv Sv Sv Sv];
            %Count the number of repetition per node
            NODES_NEL(:,ine)=NODES_NEL(:,ine)+1;
        end
    end
end
%Averaging
NODE_STRESS=NODE_STRESS./NODES_NEL;

%(c)
%Stress
s11=NODE_STRESS(1,:);
s12=NODE_STRESS(2,:);
s22=NODE_STRESS(4,:);
%stress auxiliars
tm=(s11+s22)/2;
tr=(s11-s22)/2;
tt=s12;
%vonmises
vm=tm.^2+3*(tr.^2+tt.^2);
vmM=max(vm);
figure(4)
set(4,'Position',[(nlay+1)*410 20 400 400]);
xyplot=linspace(-0.03,0.03,100);
[Xp,Yp]=meshgrid(xyplot,xyplot);
[Xr,Yr,Fr] = griddata(COORD(:,1),COORD(:,2),vm,Xp,Yp);
vlevel=[linspace(0,10*mean(vm),10)];
%contourf(Xr,Yr,Fr,vlevel)
contourf(Xr,Yr,Fr,vlevel)
colorbar
if isempty(data.hole_nodes1)
    hole=[];
else
    hole=convhull(COORD(data.hole_nodes1,:));
end
hold on
patch(COORD(data.hole_nodes1(hole),1),COORD(data.hole_nodes1(hole),2),vmM*(hole./hole)','white')
if isempty(data.hole_nodes2)
    hole=[];
else
    hole=convhull(COORD(data.hole_nodes2,:));
end
hold on
patch(COORD(data.hole_nodes2(hole),1),COORD(data.hole_nodes2(hole),2),vmM*(hole./hole)','white')
axis equal
axis tight
end
