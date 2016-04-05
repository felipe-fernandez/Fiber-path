%Read file function
% Input :
% filename  - text file obained with ABAQUS
% nlay      - number of layers
% Output:
% data      - data structure with
%  data.COORD       - Nodal coordinates
%  data.ELEM_NODE   - Table of element nodes
%  data.N_NODE      - number of nodes*degrees of freeedom
%  data.nd:         - number of nodes
%  data.N_ELEM:     - number of elements
%  data.pres:       - prescribed nodes
%  data.f_nodes1    - force nodes at hole 1
%  data.f_ele1      - force elements at hole 1
%  data.f_nodes2    - force nodes at hole 2
%  data.f_ele2      - force elements at hole 2
%  data.hole_nodes1 - nodes at the hole 1
%  data.hole_nodes2 - nodes at the hole 2
%  data.Xc          - nodal coordinate x
%  data.Yc          - nodal coordinate y

function data=recdomain(nex,ney,xL,yH)
data.nex=nex;
data.ney=ney;
data.xL=xL;
data.yH=yH;
data.h=xL/nex;

data.At=xL*yH;
%xnodes
xnodes=linspace(0,xL,nex+1);
%ynodes
ynodes=linspace(0,yH,ney+1);
%Mesh grid
[Xc,Yc]=meshgrid(xnodes,ynodes); %grid of coordinates
data.Xc=Xc;
data.Yc=Yc;
COORD=[Xc(:) Yc(:)];    %table of coordates of all the nodes
%element the correspondent nodes for each element [n1x,n1y,n2x,n2y,...,...,n4-9x,n4-9y]
ELEM_NODE=[];
for ex=1:nex
    for ey=1:ney
        %4-node bi linear element
        ELEM_NODE=[ELEM_NODE; 2*(ney+1)*(ex-1)+2*(ey-1)+[1 2 2*ney+[3 4 5 6] 3 4]];
    end
end
data.N_NODE=max(max(ELEM_NODE));  %number of nodes
data.N_ELEM=nex*ney;    %number of elements
data.pres=1:((ney+1)*2);  %prescribed degrees of freedom

%force elements in the right
data.force_ele=(ney*(nex-1)+1):(ney*nex);
data.force_node=((ney+1)*(nex)+1):((ney+1)*(nex+1));

data.nd=data.N_NODE/2;
%vectorial
data.ELEM_NODE=ELEM_NODE;

%all nodes
data.alldofs=1:data.N_NODE;
%Fixed nodes using nodeloc function for each node location
data.free=setdiff(data.alldofs,data.pres); %free degree of freedom

data.hole_nodes1=[];
data.hole_nodes2=[];

data.COORD=COORD;
data.Xc=Xc(:);
data.Yc=Yc(:);
end
