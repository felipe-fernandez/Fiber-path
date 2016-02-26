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

function data=readfile(filename)
fio=fopen(filename,'r');
linr = fgets(fio);
while linr(1)=='*'
    linr = fgets(fio);
end

ite_read=1;
while linr~='*'
    %store line
    num_lin=sscanf(linr,'%d ,%f, %f');
    %save nodal coordinates
    Xc(ite_read,1)=num_lin(2);
    Yc(ite_read,1)=num_lin(3);
    ite_read=ite_read+1;
    
    %read line
    linr=fgets(fio);
end
data.N_NODE=(ite_read-1)*2;
data.nd=data.N_NODE/2;

linr=fgets(fio);
ite_read=1;
while linr(1)~='*'
    %store line
    num_lin=sscanf(linr,'%d ,%d, %d, %d, %d');
    %save table of element nodes
    ELEMENT_NODES_s(ite_read,:)=num_lin(2:5);
    ite_read=ite_read+1;
    %read line
    linr=fgets(fio);
end
data.N_ELEM=ite_read-1;

%vectorial
ELEM_NODE=zeros(data.N_ELEM,8);
ELEM_NODE(:,1:2:7)=2*ELEMENT_NODES_s-1;
ELEM_NODE(:,2:2:8)=2*ELEMENT_NODES_s;
data.ELEM_NODE=ELEM_NODE;

linr=fgets(fio);
linr=fgets(fio);
linr=fgets(fio);
linr=fgets(fio);
linr=fgets(fio);
%prescribed degrees of freedom
data.pres=[];
pres=2*sscanf(linr,'%d ,');
data.pres=union(data.pres,union(pres,pres-1));
linr=fgets(fio);
linr=fgets(fio);
%prescribed degrees of freedom
pres=2*sscanf(linr,'%d ,');
data.pres=union(data.pres,pres);
%all nodes
data.alldofs=1:data.N_NODE;    
%Fixed nodes using nodeloc function for each node location
data.free=setdiff(data.alldofs,data.pres); %free degree of freedom


%force nodes 1(left side) 2(rigth side)
linr=fgets(fio);
linr=fgets(fio);
data.f_nodes1=[];
while linr(1)~='*'
    data.f_nodes1=union(data.f_nodes1,sscanf(linr,'%d ,'),'legacy');
    linr=fgets(fio);
end
linr=fgets(fio);
data.f_ele1=[];
while linr(1)~='*'
    data.f_ele1=union(data.f_ele1,sscanf(linr,'%d ,'),'legacy');
    linr=fgets(fio);
end
linr=fgets(fio);
data.f_nodes2=[];
while linr(1)~='*'
    data.f_nodes2=union(data.f_nodes2,sscanf(linr,'%d ,'),'legacy');
    linr=fgets(fio);
end
linr=fgets(fio);
data.f_ele2=[];
while linr(1)~='*'
    data.f_ele2=union(data.f_ele2,sscanf(linr,'%d ,'),'legacy');
    linr=fgets(fio);
end

%hole nodes
%read line
linr=fgets(fio);
data.hole_nodes1=[];
while linr(1)~='*'
    data.hole_nodes1=union(data.hole_nodes1,sscanf(linr,'%d ,'))';
    linr=fgets(fio);
end
linr=fgets(fio);
while linr(1)~='*'
    linr=fgets(fio);
end
linr=fgets(fio);
data.hole_nodes2=[];
while linr(1)~='*'
    data.hole_nodes2=union(data.hole_nodes2,sscanf(linr,'%d ,'))';
    linr=fgets(fio);
end

%coordinates
COORD=[Xc Yc];    %table of coordinates of all the nodes
data.COORD=COORD;
data.Xc=Xc;
data.Yc=Yc;
fclose(fio);

end
