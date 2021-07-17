function []=NACA4415(n,a,v)
%翼型NACA4415计算程序
%传入参数为网格分段数n(测试时建议小于100节省计算用时，实际计算建议大于100)、机翼攻角a（角度制）、来流速度v（设计参数为2）

    close all;clc;    
    %NACA4415翼型对应100m=4、10p=4、100t=15
    m=0.04;p=0.4;t=0.15;
    %设定翼型弦长c
    c=1;
    x=linspace(0,c,n);

    %建立翼型的分段中弧线方程
    syms X
    Yc_u(X)=m.*X./power(p,2).*(2*p-X./c);
    Yc_l(X)=m.*(c-X)./power(1-p,2).*(1-2*p+X./c);
    %计算中弧线在X网格上对应的点Y值
    y_c=subs(Yc_u,X,x).*(0<=x & x<p*c)+...
        subs(Yc_l,X,x).*(p*c<=x & x<=c);
    %计算翼型厚度在X网格上对应的值
    y_t=t/0.2*c.*(0.2969.*sqrt(x./c)-0.126.*(x./c)-0.3516.*power(x./c,2)+0.2843.*power(x./c,3)-0.1036.*power(x./c,4));
    %分段中弧线方程对x求导
    dYc_udx=diff(Yc_u(X));
    dYc_ldx=diff(Yc_l(X));
    %y_c求导后在X网格上对应的值
    dy_cdx=subs(dYc_udx,X,x).*(0<=x & x<p*c)+...
           subs(dYc_ldx,X,x).*(p*c<=x & x<=c);
    %中弧线上对应点夹角
    theta=atan(dy_cdx);
    %翼面划分点
    x_u=double(x-y_t.*sin(theta));
    y_u=double(y_c+y_t.*cos(theta));
    x_l=double(x+y_t.*sin(theta));
    y_l=double(y_c-y_t.*cos(theta));
    
    %上下翼面划分点合并为从前缘出发顺时针环绕的行向量
    x_surf=[x_u x_l(end-1:-1:1)];
    y_surf=[y_u y_l(end-1:-1:1)];
    
    %调用Java底层最大化窗口（如果报错无法运行整段注释掉就好）
    h = figure(1);
    jFrame = get(h,'JavaFrame');	
    set(jFrame,'Maximized',1);				
    pause(1);	
    
    %绘制翼型
    subplot(2,2,1);
    plot(x_u,y_u,x_l,y_l)
    axis equal
    title('2D NACA4415 Airfoil Profile')
    xlabel('x')
    ylabel('y')
    legend('上翼面','下翼面')

    %取翼面分段中点为控制点
    x_control_up=zeros(0,n-1);
    y_control_up=zeros(0,n-1);
    x_control_low=zeros(0,n-1);
    y_control_low=zeros(0,n-1);
   
    for i=1:n-1
        %上翼面控制点
        x_control_up(i)=(x_u(i)+x_u(i+1))/2;
        y_control_up(i)=(y_u(i)+y_u(i+1))/2;
        %下翼面控制点
        x_control_low(i)=(x_l(i)+x_l(i+1))/2;
        y_control_low(i)=(y_l(i)+y_l(i+1))/2;
    end
    %上下翼面控制点合并为从前缘出发顺时针环绕的行向量
    x_ctrl=[x_control_up x_control_low(end:-1:1)];
    y_ctrl=[y_control_up y_control_low(end:-1:1)];
    for i=1:2*(n-1)
        l(i)=sqrt(x_ctrl(i)^2+y_ctrl(i)^2);
    end
    
    length=zeros(1,2*(n-1));
    for i=1:2*(n-1)
        % 每段翼面微元的长度
        length(i)=sqrt((x_surf(i+1)-x_surf(i))^2+(y_surf(i+1)-y_surf(i))^2);
    end
    
    [r,~,~]=rSolver(n,a,v,x_surf,y_surf,x_ctrl,y_ctrl,length,l);
    
    %切向速度行向量
    u=(r(1:2*(n-1),1))';
    %求压力系数行向量
    Cp=1-(u./v).^2;
    
    %压力系数曲线
    subplot(2,2,2);
    plot(x_control_up,Cp(1:n-1),'.',x_control_low(1:(n-2)),Cp(2*n-2-1:-1:n),'.')
    title('Pressure Coefficient')
    xlabel('x')
    ylabel('压力系数Cp')
    legend('上翼面','下翼面')
    set(gca,'YDir','reverse'); 
    
    %升力系数曲线
    subplot(2,2,3);
    cl=[];
    cmle=[];
    bb=[-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10];

    for q=-5:10
       [~,T1,T2]=rSolver(n,q,v,x_surf,y_surf,x_ctrl,y_ctrl,length,l);
       Cl=2*T1/(v*c);
       Cmle=2*T2/(v*c);
       cl(q+6)=Cl;
       cmle(q+6)=Cmle;
    end

    %薄翼理论
    syms z;
    q=int((0.25*cos(z)-0.05)*(cos(z)-1),0,acos(0.2));
    b=int((5/45*cos(z)-1/45)*(cos(z)-1),acos(0.2),pi);
    q0=-1/pi*(q+b);
    e=int((0.25*cos(z)-0.05)*(cos(2*z)-cos(z)),0,acos(0.2));
    f=int((5/45*cos(z)-1/45)*(cos(2*z)-cos(z)),acos(0.2),pi);
    cml0=(1/2)*(e+f);
    z=-5:1:10;
    cmle1=-(2*pi*(z/180*pi-q0))/4+cml0;
    cl1=2*pi*(z/180*pi-q0);

    plot(z,cl1)
    hold on;
    plot(bb,cl); 
    hold off
    title("升力系数曲线")  
    legend('薄翼理论','环量理论')

    %力矩系数曲线
    subplot(2,2,4);
    plot(z,cmle1);
    hold on 
    plot(bb,cmle);
    title("阻力系数曲线")  
    hold off
    legend('薄翼理论','环量理论')

     %作出时间线
    flowSim(n,a,v,u,x_ctrl,y_ctrl,x_surf,y_surf,length)
end
