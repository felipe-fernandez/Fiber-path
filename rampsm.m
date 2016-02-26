%step function
function [h,dh]=rampsm(phi,p)
% %power law ramp smooth function p-> infty then h is sharp ramp
% h=((((phi+1)).^p+1).^(1/p))-1;
% h(phi<-1)=0;
% %derivative of the ramp function
% dh=((phi + 1).^(p - 1)).*((phi + 1).^p + 1).^(1/p - 1);
% dh(phi<-1)=0;

%Piece wise ramp function
a=10^(-p);
b=10^(0);
h=a*phi+b*(phi-1).^2;
h(phi<1)=a*phi(phi<1);
h(phi<0)=0;

dh=a+2*b*(phi-1);
dh(phi<1)=a;
dh(phi<0)=0;
end
