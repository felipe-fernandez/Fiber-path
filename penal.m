%Penalization of the density
function [rho,drho]=penal(chi)
p=1;
%penalization
rho=chi^p;
rho(rho>1)=1;
rho(rho<0)=0;
if chi>1 || chi<0
    drho=0;
else
    drho=p*chi^(p-1);
end

% rho=1;
% drho=0;

% epsilon=.5;
% chit=chi-.5;
% %heaviside function
% rho=((3/4)*(chit/epsilon-(chit.^3)/(3*(epsilon^3)))+1/2);
% rho(chit>epsilon)=1;
% rho(chit<-epsilon)=0;
% %derivative
% drho=3/4*(1/epsilon-(3*chit.^2)/(3*epsilon^3));
% drho(chit>epsilon)=0;
% drho(chit<-epsilon)=0;

end
