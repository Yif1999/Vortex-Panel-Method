function [r,T1,T2]=rSolver(n,a,v,x_surf,y_surf,x_ctrl,y_ctrl,length,l)
%涡强求解器
%传参分别为网格分段数n、机翼攻角a、来流速度v、一面信息
    
    %攻角调整
    [x_surf,y_surf]=angleChange(x_surf,y_surf,a);
    [x_ctrl,y_ctrl]=angleChange(x_ctrl,y_ctrl,a);	

    %创建方程组矩阵
    A=zeros(2*n-1,2*n-2+1);
    B=zeros(2*n-1,1);
    for i=1:2*(n-1)
        B(i)=v*y_ctrl(i);
        for j=1:2*(n-1)           
            A(i,j)=Calc(x_ctrl(i),y_ctrl(i),x_surf(j),y_surf(j),x_surf(j+1),y_surf(j+1));
        end
        A(i,2*(n-1)+1)=1;
    end
    A(2*(n-1)+1,n-1:n)=1;
    A(2*(n-1)+1,2*(n-1)+1)=0;

    %求解面源涡强
    r=(-1)*A\B;
    
    %环量计算
    for i=1:2*(n-1)
    T1(i)=r(i)*length(i);
    end
    T1=sum(T1);
    
    %俯仰力矩系数计算
    for i=1:2*(n-1)
    T2(i)=l(i)*r(i)*length(i);
    end
    T2=-sum(T2);
    
 end