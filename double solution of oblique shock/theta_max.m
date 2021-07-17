function [thetaM] = theta_max(M1,betaM)
%THETAM ¼ÆËãtheta_max
    betaM=betaM/180*pi;
    numerator=2*((M1^2-1)*tan(betaM)^2-1);
    denominator=tan(betaM)*((1.4*M1^2+2)*(1+tan(betaM)^2)+M1^2*(1-tan(betaM)^2));
    thetaM=atan(numerator/denominator)/pi*180;
end

