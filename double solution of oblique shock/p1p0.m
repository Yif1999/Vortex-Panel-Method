function [p10]=p1p0(M0,beta1)
%根据Ma0和beta1求p1/p0
    gama=1.4;
    p10=2*gama*M0^2*(sin(beta1*pi/180))^2/(gama+1)-(gama-1)/(gama+1);
end