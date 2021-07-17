function [thetaD] = thetaD(theta1)
%THETAD µü´úÇó½âthetaD
    M=4;
    theta2=-1;
    while abs(theta2-theta1)>0.001
        if theta2~=-1
            theta1=(theta1+theta2)/2;
        end
        beta1=beta(M,theta1);
        M1=Ma(M,beta1);
        beta2=beta_max(M1);
        theta2=theta(M1,beta2);
    end
    thetaD=theta2;
end