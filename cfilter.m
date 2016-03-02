% Cone filter function

function [G,data]=cfilter(data,efilter)

%get data
nlay=data.nlay;
ELEM_NODE=data.ELEM_NODE;
COORD=data.COORD;

nodev=1:(data.nd);
elev=1:(data.N_ELEM);
elem_nod=ELEM_NODE(:,2:2:end)/2;
nd=data.nd;
ne=data.N_ELEM;

Ngauss=1;   %Center of each element
[egv,wg] = GLTable(Ngauss);
Area_el=zeros(ne,1);
xc_el=zeros(ne,2);

%All the elements
for ele=elev
    %Index of element nodes scalar field
    iens=elem_nod(ele,:);

    %nodal coordinates of the element
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
            %center
            xc_el(ele,:)=Ne*Xe';
            %area of each element
            Area_el(ele)=Area_el(ele)+wg(eit)*wg(nit)*Jdet;
        end
    end
end
%store for volume constraint
data.Area_el=Area_el;
%store center of the elements for plot
data.xc_el=xc_el;

if efilter==0
    Gc=eye(nd,nd);
else
    %initialize parameters
    nd=data.nd;
    Gc=sparse(nd,nd);
    %loop over elements
    for nodp=nodev
        %nodal coordinate
        xo=COORD(nodp,:);
        %distance from the node to each other node
        disn=sqrt((COORD(:,1)-xo(1)).^2+(COORD(:,2)-xo(2)).^2);
        %nodes
        nodesum=nodev(disn<efilter);
        %elements
        elesum=union(union(elev(ismember(elem_nod(:,1),nodesum)),...
            elev(ismember(elem_nod(:,2),nodesum))),...
            union(elev(ismember(elem_nod(:,3),nodesum)),...
            elev(ismember(elem_nod(:,4),nodesum))));
        %center of element in vector
        xe=xc_el(elesum,:);
        %distance
        dise=sqrt((xe(:,1)-xo(1)).^2+(xe(:,2)-xo(2)).^2);
        %index of nodes
        iene=elem_nod(elesum,:);
        %filter
        filtera=Area_el(elesum).*(efilter-dise)/(efilter);
        filtera(dise>efilter)=0;
        %filter
        Gc(nodp,iene(:,1))=Gc(nodp,iene(:,1))+filtera';
        Gc(nodp,iene(:,2))=Gc(nodp,iene(:,2))+filtera';
        Gc(nodp,iene(:,3))=Gc(nodp,iene(:,3))+filtera';
        Gc(nodp,iene(:,4))=Gc(nodp,iene(:,4))+filtera';
        
    end
    
    %average
    for nod=1:data.nd
        %average
        Gc(nod,:)=Gc(nod,:)/sum(Gc(nod,:));
    end
end

%layers filter
G=sparse((nd+ne)*nlay,(nd+ne)*nlay);
for lay=1:nlay
    G((lay-1)*nd+(1:nd),(lay-1)*nd+(1:nd))=Gc;
end
%filter for elements
Ge=eye(nlay*ne);
G((nd*nlay+1):end,(nd*nlay+1):end)=Ge;
end