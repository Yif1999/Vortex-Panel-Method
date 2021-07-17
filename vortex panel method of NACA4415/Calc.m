function sum=Calc(x_c,y_c,x1,y1,x2,y2)
    n=300;
    x=linspace(x1,x2,n);
    dx=(x2-x1)/n;
    dy=(y2-y1)/n;
    dl=sqrt((x1-x2)^2+(y1-y2)^2)/n;
    sum=0;
    for i=1:n
        xi=x1+i*dx-dx/2;
        yi=y1+i*dy-dy/2;
        vector_d=x_c-xi+1i*(y_c-yi);
        c=log(abs(vector_d));
        temp=dl/(2*pi)*c;
        sum=sum+temp;
    end
end