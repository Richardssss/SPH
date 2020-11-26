clear;
%测试10*20板碰撞
%单位m，s
%建立模型：将平板分为10*20网格，每个网格内有一个粒子，进行离散处理
%左
xl = linspace(-0.02,-0.01,11);
yl = linspace(-0.01,0.01,21);
%右
xr = linspace(0,0.01,11);
yr = linspace(-0.01,0.01,21);

%2*10*20个粒子的坐标
x = zeros(20,20);
y = zeros(20,20);
%左侧平板，动板
for i = 1:20
    for j =1:10
        x(i,j) = ( xl(j) + xl(j+1) ) / 2;
        y(i,j) = ( yl(i) + yl(i+1) ) / 2;
    end
end
%右侧静板
for i = 1:20
    for j = 11:20
        x(i,j) = ( xr(j-10) + xr(j-10+1) ) / 2;
        y(i,j) = ( yr(i) + yr(i+1) ) / 2;
    end
end

%粒子的速度
vx = zeros(20,20);
vy = zeros(20,20);
%左侧板500m/s
for i = 1:20
    for j =1:10
        vx(i,j) = 500;
    end
end
%初始化其他参数
rho0 = 18.087e+3;
m = rho0 * 0.001^2;
rho = zeros(20,20);
for i = 1:20
    for j = 1:20
        rho(i,j) = rho0;
    end
end

e = zeros(20,20);
for i = 1:20
    for j = 1:20
        e(i,j) = 0;
    end
end

c0 = 2510;
s = 1.51;
Gamma = 2.13;

PH = (c0^2*(1/rho0 - 1./rho))./(1/rho0 - s*(1/rho0 - 1./rho)).^2;

eH = 0.5*PH.*(1/rho0 - 1./rho);

P = PH + Gamma.*rho.*(e - eH);
% sigma = -P;
figure;
subplot(2,2,1);
plot(x,y,'o');
xlabel('x');
ylabel('y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 22;
dt = 1e-6;
h = 0.001*2;
for t = 1:T
    P1 = P;
    x1 = x;
    y1 = y;
    rho1 = rho;
    vx1 = vx;
    vy1 = vy;
    %先算x的变化
    x = x + vx*dt;
    y = y + vy*dt;
    for iy = 1:20
        for ix = 1:20
            drho = 0;
            de = 0;
            dvx = 0;
            dvy = 0;
            for jy = 1:20
                for jx = 1:20
                    Rx = x1(iy,ix)-x1(jy,jx);
                    Ry = y1(iy,ix)-y1(jy,jx);
                    R = sqrt(Rx^2+Ry^2)/h;
                    if R<=1 && R>0
                        [dWx,dWy] = dW(Rx,Ry,h);
                        drho = drho + m*(vx1(iy,ix)-vx1(jy,jx))*dWx + m*(vy1(iy,ix)-vy1(jy,jx))*dWy;
                        de = de + 0.5*(-P1(iy,ix)./rho1(iy,ix)^2-P1(jy,jx)./rho1(jy,jx)^2)*(m*(vx1(iy,ix)-vx1(jy,jx))*dWx + m*(vy1(iy,ix)-vy1(jy,jx))*dWy);
                        dvx = dvx + (-P1(iy,ix)./rho1(iy,ix)^2-P1(jy,jx)./rho1(jy,jx)^2)*m*dWx;
                        dvy = dvy + (-P1(iy,ix)./rho1(iy,ix)^2-P1(jy,jx)./rho1(jy,jx)^2)*m*dWy;
                    end
                end
            end
            rho(iy,ix) = rho(iy,ix) + drho*dt;
            e(iy,ix) = e(iy,ix) + de*dt;
            vx(iy,ix) = vx(iy,ix) + dvx*dt;
            vy(iy,ix) = vy(iy,ix) + dvy*dt;
        end
    end
    %求e
   
    %求PH;c0 = 0.251;s = 1.51;
    PH = (c0^2*(1/rho0 - 1./rho))./(1/rho0 - s*(1/rho0 - 1./rho)).^2;
    %求eH
    eH = 0.5*PH.*(1/rho0 - 1./rho);
    %求P
    P = PH + Gamma.*rho.*(e - eH);
end
subplot(2,2,2);
plot(x,y,'o');
xlabel('x');
ylabel('y');
subplot(2,2,3);
plot(x(10,:),P(10,:));
xlabel('x');
ylabel('P(-sigma)');
subplot(2,2,4);
plot(x(10,:),vx(10,:));
xlabel('x');
ylabel('vx');