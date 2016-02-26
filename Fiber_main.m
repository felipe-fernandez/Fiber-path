% Additive manufacturing tool path optimization code
% Author: Felipe Fernandez
% Main function that do all the magic
function Fiber_main()

%cleaning the house first
clear all
clc
close all
%Ready to get nice plots
% Change default axes fonts.
%set(0,'DefaultAxesFontName', 'Helvetica')
set(0,'DefaultAxesFontSize', 10)
% Change default text fonts.
%set(0,'DefaultTextFontname', 'Helvetica')
set(0,'DefaultTextFontSize', 10)
set(0, 'defaultTextInterpreter', 'latex');
set(0,'defaultaxesfontsize',10)
%set(0,'defaultaxesfontname','Helvetica','defaultaxesfontsize',10)


% Parameters and constants
th=5/1000;      %Thickness (m)
loadv=0;        %traction force magnitude for applied nodes
fpres=[1e6;0];  %prescribed pressure force for applied nodes
nmin=.1;        %minimun gradient constraint
nmax=10;        %maximum gradient constraint
rmin=1e-3;      %minimun radius for toolpath
pow=2;          %power that smooth ramp function
efilter=5e-3;   %radius of filter
cdiv=10;        %divergence constraint parameter
rc='g';         %requested constraints 'g' gradient 'k' curvature 'd' divergence 
nlay=2;         %number of layers

% Read mesh file
filename='fex2_2.txt';

% Initialization
[data,UG0,FG,G]=initializationf(filename,nlay,loadv,th,rc,efilter,fpres);

% %Topology optimization procedure
% %Initial guess
% %dv=0.9*(data.Yc+15/1000);%rand(size(data.Yb));%data.Yb,yH-data.Yb);%min(min(data.Yb,data.Xb),min(yH-data.Yb,xL-data.Xb));
% dv=sqrt(data.Yc.^2+data.Xc.^2);
% %initial fea to normalize
% [theta,dtheta]=feafun(dv,G,data,UG0,FG,th,1,nameplot,nmax);
% [cons,ceq,dcons,dceq]=nlcn(dv,G,data,nmin,nmax,pow,rmin);
% 
% cons
% %finite differences
% ddiff=1e-6;
% for var=1:10
%     alphad=dv;
%     alphad(var)=alphad(var)+ddiff;
%     %finite element and post processing function
%     [thetad,dthetad]=feafun(alphad,G,data,ELEM_NODE,UG0,FG,COORD,th,1,nameplot,nmax);
%     [consd,ceqd,dconsd,dceqd]=nlcn(alphad,G,data,nmin,nmax,pow,rmin);
%     %derivative
%     dtheta_fd(var)=(thetad-theta)/ddiff;
%     dcons_fd(:,var)=(consd-cons)/ddiff;
% end
% dtheta_fd
% dtheta(1:10)
% dcons_fd
% dcons(1:10,:)'
% afdasdfadsf

%initial level set function
dv0=data.Yc;
dv=[];
for i=1:nlay
    dv=[dv;dv0]/0.03;
end

%initial fea to normalize
[c0,dtheta]=feafun(dv,eye(nlay*data.nd),data,UG0,FG,th,1,nmax);
[cons,ceq,dcons,dceq]=nlcn(dv,eye(nlay*data.nd),data,nmin,nmax,pow,rmin,cdiv,rc);

% dv=(-sqrt(19/20)*data.Yc+sqrt(1/20)*data.Xc);
% if nlay>1
%     dv=[dv;(sqrt(19/20)*data.Yc+sqrt(1/20)*data.Xc)];
% end
% if nlay==3
%     dv=[dv;(1+.2*rand(data.nd,1)).*data.Yc];
% end

dv=data.Yc;
if nlay>1
    dv=[dv;data.Xc];
end
if nlay==3
    dv=[dv;(1+.2*rand(data.nd,1)).*data.Yc];
end
dv=dv/0.03;


%finite element and post processing function
obFUN=@(dalpha)feafun(dalpha,G,data,UG0,FG,th,c0,nmax);
conFUN=@(dalpha)nlcn(dalpha,G,data,nmin,nmax,pow,rmin,cdiv,rc);
%Gradients of objective and constraints are given by the user
%algorithm is set as interior-point as suggest matlab help
options = optimset('GradObj','on',...
    'Display','on','GradConstr','on',...
    'Tolfun',1e-3,'TolCon',1e-3,'TolX',1e-4,'MaxIter',3,...
    'Algorithm','interior-point','Display','iter','MaxFunEvals',500,'AlwaysHonorConstraints','bounds');
%call fmincon function embedded in Matlab
[dvo,fval] =fmincon(obFUN,dv,[],[],[],[],-ones(nlay*data.nd,1),ones(nlay*data.nd,1),conFUN,options);

[theta,dtheta]=feafun(dvo,G,data,UG0,FG,th,c0,nmax);
[cons,ceq,dcons,dceq]=nlcn(dvo,G,data,nmin,nmax,pow,rmin,cdiv,rc);
disp(['c=' num2str(theta,'%1.8f') ])
save([data.nameplot '.mat']);
end


