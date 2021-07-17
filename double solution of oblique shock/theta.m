function [theta] = theta(M,beta)
%THETA ¸ù¾İbetaËãtheta
    beta=2*pi/360*beta;
    theta=atan(2*cot(beta)*(M^2*sin(beta)^2-1)/(M^2*(1.4+cos(2*beta))+2))/pi*180;
end

