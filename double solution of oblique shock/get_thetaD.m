function [thetaN]=get_thetaD(M0)
%»ñµÃthetaD
    beta0=30;
    theta0=theta(M0,beta0);
    while(1)
        M1=Ma(M0,beta0);
        beta1=beta_max(M1);
        theta1=theta_max(M1,beta1);
        if(abs(theta1-theta0)<=0.02)
            break;
        else
            theta0=(theta1+theta0)/2;
            beta0=beta(M0,theta0);
        end
    end
    thetaN=theta0;
end