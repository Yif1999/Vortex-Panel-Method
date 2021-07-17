function thetaN=thetaN(M0)
%»ñµÃthetaN
    beta0=30;
    theta0=theta(M0,beta0);
    while(1)
        M1=Ma(M0,beta0);
        p30=p3p0(M0);
        p10=p1p0(M0,beta0);
        beta1=getBetaFromP(M1,p10,p30);
        theta1=theta(M1,beta1);
        if(abs(theta1-theta0)<=0.02)
            break;
        else
            theta0=(theta1+theta0)/2;
            beta0=beta(M0,theta0);
    end
    thetaN=theta0;
end