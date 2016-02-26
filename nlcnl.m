%non linear constraints for each layer
function [consr,dconsr,vfrac]=nlcnl(dv,data,ELEM_NODE,nmin,nmax,pow,rmin,lay,nlay,COORD,cdiv,rc)
cons=[0;0;0;0];
dcons=zeros(4,nlay*data.nd);

Ngauss=1;   %Number of Gauss points
[egv,wg] = GLTable(Ngauss);

%initialize values
area=0;
vfrac=0;
dvfrac=zeros(1,nlay*data.nd);
DPHIv=zeros(data.nd,2);
dDPHIvx=sparse(data.nd,data.nd);
dDPHIvy=sparse(data.nd,data.nd);
counv=zeros(data.nd,1);
nd=data.nd;

%integration in the domain
for ele=1:data.N_ELEM
    %Index of element nodes vector field
    ienv=ELEM_NODE(ele,:);
    %index of nodes scalar field
    ien=ienv(2:2:end)/2;
    
    Xe=COORD(ien,:)';
    %Gauss quadrature integration
    for eit=1:Ngauss
        eg=egv(eit);
        for nit=1:Ngauss
            ng=egv(nit);
            
            %Derivative of the shape function
            dNB=1/4*[-(1-ng) 1-ng 1+ng -(1+ng);...
                -(1-eg) -(1+eg) 1+eg 1-eg]';
            J=Xe*dNB;
            
            %Derivative of the shape function
            Be=(dNB/J)';
            
            %gradient level set function
            dphi=Be*dv(ien);
            
            %normal of the gradient (diferentiable at 0)
            Normdphi=sqrt((dphi'*dphi));
            %derivative of the gradient
            dNormdphi=dphi'/sqrt(dphi'*dphi)*Be;
            
            
            %extrapolate to the nodes
            DPHIv(ien,1)=DPHIv(ien,1)+dphi(1);
            DPHIv(ien,2)=DPHIv(ien,2)+dphi(2);
            
            dDPHIvx(ien,ien)=dDPHIvx(ien,ien)+...
                [Be(1,:); Be(1,:); Be(1,:); Be(1,:)];
            dDPHIvy(ien,ien)=dDPHIvy(ien,ien)+...
                [Be(2,:); Be(2,:); Be(2,:); Be(2,:)];
            counv(ien)=counv(ien)+1;
        end
    end
end

Mcoun=repmat(counv,1,nd);
%Gradient 
DPHIv=DPHIv./[counv counv];
dDPHIvx=dDPHIvx./Mcoun;
dDPHIvy=dDPHIvy./Mcoun;

%normal vector
nDPHIv=sqrt(DPHIv(:,1).^2+DPHIv(:,2).^2);
norv=DPHIv./([nDPHIv nDPHIv]);

dnorvx=sparse(data.nd,data.nd);
dnorvy=sparse(data.nd,data.nd);
for jnod=1:data.nd
    dnDPHIv=(DPHIv(jnod,1)*dDPHIvx(jnod,:)+DPHIv(jnod,2)*dDPHIvy(jnod,:))/nDPHIv(jnod);
    dnorvx(jnod,:)=dnorvx(jnod,:)+dDPHIvx(jnod,:)/nDPHIv(jnod)-DPHIv(jnod,1)*dnDPHIv/(nDPHIv(jnod)^2);
    dnorvy(jnod,:)=dnorvy(jnod,:)+dDPHIvy(jnod,:)/nDPHIv(jnod)-DPHIv(jnod,2)*dnDPHIv/(nDPHIv(jnod)^2);
end


Ngauss=3;   %Number of Gauss points
[egv,wg] = GLTable(Ngauss);
NODE_val=zeros(1,nd);
NODE_nel=zeros(1,nd);

%integration in the domain
for ele=1:data.N_ELEM
    %Index of element nodes vector field
    ienv=ELEM_NODE(ele,:);
    %index of nodes scalar field
    ien=ienv(2:2:end)/2;
    %index global desgin variable
    ienl=ien+(lay-1)*nd;
    Xe=COORD(ien,:);
    %Gauss quadrature integration
    for eit=1:Ngauss
        eg=egv(eit);
        for nit=1:Ngauss
            ng=egv(nit);
            
            %Derivative of the shape function
            dNB=1/4*[-(1-ng) 1-ng 1+ng -(1+ng);...
                -(1-eg) -(1+eg) 1+eg 1-eg];
            ddNB=[0 0 0 0;1 -1 1 -1;0 0 0 0];
            %jacobian
            J=dNB*Xe;
            Jdet=det(J);
            %Derivative of the shape function
            Be=J\dNB;
            %Le [d2x/de2 d2y/de2;d2x/dedn d2y/dedn;d2x/dn2 d2y/dn2]
            Le=ddNB*Xe;
            
            %Se [(dx/de)^2 dxdx/dedn (dx/dn)^2;
            %2dxdy/dede dxdy/dedn+dxdy/dnde 2dxdy/dndn;
            %(dy/de)^2 dydy/dedn (dy/dn)^2]^T
            Se=[J(1,1)^2 J(1,1)*J(2,1) J(2,1)^2;...
                2*J(1,1)*J(1,2) J(1,1)*J(2,2)+J(1,2)*J(2,1) 2*J(2,1)*J(2,2);...
                J(1,2)^2 J(1,2)*J(2,2) J(2,2)^2]';
            %Second derivative Shape functions
            Ee=Se\(ddNB-Le*Be);
                                   
            %gradient level set function
            dphi=Be*dv(ien);
            d2phi=Ee*dv(ien);
            
            %function
            c4=(d2phi(1)+d2phi(3))^2;
            dc4=2*(d2phi(1)+d2phi(3))*(Ee(1,:)+Ee(3,:));
            
            divphi=Be(1,:)*DPHIv(ien,1)+Be(2,:)*DPHIv(ien,2);
            c5=(divphi)^2;
            dc5=2*(divphi)*(Be(1,:)*dDPHIvx(ien,:)+Be(2,:)*dDPHIvy(ien,:));
            
            %curvature as gradient of the normal vector
            kur=Be(1,:)*norv(ien,1)+Be(2,:)*norv(ien,2);
            dkur=Be(1,:)*dnorvx(ien,:)+Be(2,:)*dnorvy(ien,:);
            
            %normal of the gradient (diferentiable at 0)
            Normdphi=sqrt((dphi'*dphi));
            %derivative of the gradient
            dNormdphi=dphi'/sqrt(dphi'*dphi)*Be;
            
            %norm of gradient square
            Normdphi2=Normdphi^2;
            %derivative of the normal of gradient
            dNormphi2=2*Normdphi*dNormdphi;
            
            
            
            %constraint fibers shouldnt be far
            [c1,dc1]=rampsm(2-Normdphi2/(nmin^2),pow);
            dc1=dc1*dNormphi2/(nmin^2);
            %constraint fibers shouldnt be too close (overlap)
            [c2,dc2]=rampsm(Normdphi2/(nmax^2),pow);
            dc2=dc2*dNormphi2/(nmax^2);
            %constraint in curvature
            [c3,dc3]=rampsm((kur*rmin)^2,pow);
            dc3=dc3*2*kur*(rmin^2)*dkur;
            
            if abs(kur*rmin)>=1
               %disp(['r is smaller than rmin, r is ' num2str(1/kur)])
            end
            %area
            area=area+wg(eit)*wg(nit)*Jdet;
            
            %function
            cons=cons+wg(eit)*wg(nit)*[c1;c2;c3;c5]*Jdet;
            %dcons(1:3,ienl)=dcons(1:3,ienl)+wg(eit)*wg(nit)*[-dc1;dc2;dc4]*Jdet;
            dcons(1:2,ienl)=dcons(1:2,ienl)+wg(eit)*wg(nit)*[-dc1;dc2]*Jdet;
            dcons(3:4,(1:nd)+(lay-1)*nd)=dcons(3:4,(1:nd)+(lay-1)*nd)+wg(eit)*wg(nit)*[dc3;dc5]*Jdet;
            %volume fraction
            chi=Normdphi/nmax;
            
            %nodal averaging
            NODE_val(:,ien)=...
                NODE_val(:,ien)+[c5];
            %Count the number of repetition per node
            NODE_nel(:,ien)=NODE_nel(:,ien)+1;
        end
    end
end

%ramp at phi=0
h0=rampsm(1,pow);
% %normalized by the area of ramp at phi=0
% cons(1:2)=(cons(1:2)/(data.At*h0))-1;
% cons(3)=cons(3)/(data.At*h1)-1;
% dcons=[dcons(1:2,:)'/(data.At*h0) dcons(3,:)'/(data.At*h1)];
% %cons=(cons(1:2)/(data.At*h0))-1;
% %dcons=[dcons(1:2,:)'/(data.At*h0)];


% cons=cons(3);
% dcons=dcons(3,:)';

cons(1:3)=(cons(1:3)/(data.At*h0))-1;
cons(4)=cons(4)/(data.At*cdiv^2)-1;
dcons=[dcons(1:3,:)'/(data.At*h0) dcons(4,:)'/(data.At*cdiv^2)];

consr=cons(ismember('ggkd',rc));
dconsr=dcons(:,ismember('ggkd',rc));


% %Averaging
% NODE_val=NODE_val./NODE_nel;
% 
% %(c)
% val=NODE_val;
% vmM=max(val);
% figure(lay+4)
% set(lay+4,'Position',[10+(lay-1)*410 420 400 400]);
% xyplot=linspace(-0.03,0.03,100);
% [Xp,Yp]=meshgrid(xyplot,xyplot);
% [Xr,Yr,Fr] = griddata(COORD(:,1),COORD(:,2),val,Xp,Yp);
% %contourf(Xr,Yr,Fr,vlevel)
% contourf(Xr,Yr,Fr)
% colorbar
% if isempty(data.hole_nodes1)
%     hole=[];
% else
%     hole=convhull(COORD(data.hole_nodes1,:));
% end
% hold on
% patch(COORD(data.hole_nodes1(hole),1),COORD(data.hole_nodes1(hole),2),vmM*(hole./hole)','white')
% if isempty(data.hole_nodes2)
%     hole=[];
% else
%     hole=convhull(COORD(data.hole_nodes2,:));
% end
% hold on
% patch(COORD(data.hole_nodes2(hole),1),COORD(data.hole_nodes2(hole),2),vmM*(hole./hole)','white')
% axis equal
% axis tight
% hold off
end
