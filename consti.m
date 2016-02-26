%Constitutive relationship in matrix form
function matC=consti()
% EL=7.689e9;
% ET=3.4e9;
% GLT=2.695e9;
% nuTL=0.13709;
% nuLT=nuTL*EL/ET

EL=137.9e9;
ET=10.34e9;
nuLT=0.29;
GLT=6.89e9;
nuTL=ET*nuLT/EL;
%orthotropic material in plane stress
c1111=EL/(1-nuLT*nuTL);
c2222=EL*nuTL/(1-nuLT*nuTL);
c1122=ET/(1-nuLT*nuTL);
c2121=GLT;

%constittive matrix
matC=[c1111 0 0 c1122;...
    0 c2121 c2121 0;...
    0 c2121 c2121 0;...
    c1122 0 0 c2222];

% E=80e6;
% nu=0.3;
% %Lame constants
% lam=E*nu/(1+nu)/(1-2*nu);   %plane strain
% miu=E/(2*(1+nu));
% lam=2*miu*lam/(lam+2*miu);  %plane stress
% %Constitutive matrix isotropic
% matC=[lam+2*miu 0 0 lam;0 miu miu 0;0 miu miu 0;lam 0 0 lam+2*miu];

end
