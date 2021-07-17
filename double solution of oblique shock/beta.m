function [beta] = beta(M,theta)
% 给出M和theta求beta
    theta=theta/180*pi;
    beta=fsolve(@(B)myfun(M,B,theta),theta);
    beta=beta/pi*180;
end

function F = myfun(M,B,theta)
F = 2*cot(B)*(M^2*sin(B)^2-1)/(M^2*(1.4+cos(2*B))+2)-tan(theta);
end
