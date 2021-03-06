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
nmin=1e-10;        %minimun gradient constraint
nmax=10;        %maximum gradient constraint
rmin=2e-3;      %minimun radius for toolpath
pow=2;          %power that smooth ramp function
efilter=2.5e-3;   %radius of filter
cdiv=1;        %divergence constraint parameter
rc='gk';         %requested constraints 'g' gradient 'k' curvature 'd' divergence 
nlay=2;         %number of layers
vfc=.5;        %volume fraction constraint
vel='c';        %velocity options: 'c' constant 'v' variable
numv=5;         %number of random variables to check finite differences 
prob='c';       %problem to solve: 'r' rectangle, 'c' clevis

% Read mesh file
filename='fex2_2.txt';

%just plot result from a .mat file
%plotfile(nlay, rc, efilter, vel, prob)
%return 

% Initialization
[data,UG0,FG,G]=initializationf(filename,nlay,th,rc,efilter,vel,prob);

%initial level set function
dv0=vfc*(data.Yc)*nmax;
%dv1=sqrt(data.Yc.^2+data.Xc.^2);
dv=repmat(dv0,nlay,1)/.03;
%dv=[dv0;dv1]/.03;
%initial element densities
if vel~='c'
    dv=[dv;vfc*ones(data.N_ELEM*nlay,1)];
end

%initial fea to normalize
[c0,dtheta]=feafun(dv,eye(nlay*(data.nd+data.N_ELEM)),data,UG0,FG,th,1,nmax,1,vel,prob);
[cons,ceq,dcons,dceq]=nlcn(dv,eye(nlay*(data.nd+data.N_ELEM)),data,nmin,nmax,pow,rmin,cdiv,rc,vfc,vel);

fac=sqrt(2)/2;
dv=fac*data.Yc+fac*data.Xc;

if nlay>1
    dv=[dv;fac*data.Yc-fac*data.Xc];
end
if nlay==3
    dv=[dv;(1+.2*rand(data.nd,1)).*data.Yc];
end
dv=vfc*nmax*dv/0.03;
%initial element densities
if vel~='c'
    dv=[dv;vfc*ones(data.N_ELEM*nlay,1)];
end
%finite element and post processing function
obFUN=@(dalpha)feafun(dalpha,G,data,UG0,FG,th,c0,nmax,1,vel,prob);
conFUN=@(dalpha)nlcn(dalpha,G,data,nmin,nmax,pow,rmin,cdiv,rc,vfc,vel);

% %check sensitivities with finite differences
% checkfindiff(obFUN,conFUN,data,numv,vel);
% return

%Gradients of objective and constraints are given by the user
%algorithm is set as interior-point as suggest matlab help
options = optimset('GradObj','on',...
    'Display','on','GradConstr','on',...
    'Tolfun',1e-4,'TolCon',1e-4,'TolX',1e-4,'MaxIter',1000,...
    'Algorithm','interior-point','Display','iter','MaxFunEvals',1000,'AlwaysHonorConstraints','bounds');
%call fmincon function embedded in Matlab
low=-20*ones(nlay*data.nd,1);
if vel~='c'
    low=[low;1e-3*ones(nlay*data.N_ELEM,1)];
end
[dvo,fval] =fmincon(obFUN,dv,[],[],[],[],low,20*ones(nlay*(data.nd+(vel~='c')*data.N_ELEM),1),conFUN,options);

[theta,dtheta]=feafun(dvo,G,data,UG0,FG,th,c0,nmax,1,vel,prob);
[cons,ceq,dcons,dceq]=nlcn(dvo,G,data,nmin,nmax,pow,rmin,cdiv,rc,vfc,vel);
disp([data.nameplot '  g=' num2str(nmin) '/' num2str(nmax) '  filter=' num2str(efilter) '  rmin=' num2str(rmin)  '  cdiv='  num2str(cdiv) '  c=' num2str(theta,'%1.8f') ])
save([data.nameplot '.mat']);
end


