function [betaM]=beta_max(M1)
%求解最大beta
    gama=1.4;
    sinsq=(((gama+1)/4)*M1^2-1+((1+gama)*(1+((gama-1)*M1^2)/2+((gama+1)*M1^4)/16))^0.5)/(gama*M1^2);
    betaM=asin(sinsq^0.5)*180/pi;
end