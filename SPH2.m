clear;
%100*200�����ײ
%��λm��s
%����ģ�ͣ���ƽ���Ϊ100*200����ÿ����������һ�����ӣ�������ɢ����
%��
xl = linspace(-0.02,-0.01,101);
yl = linspace(-0.01,0.01,201);
%��
xr = linspace(0,0.01,101);
yr = linspace(-0.01,0.01,201);

%2*100*200�����ӵ�����
x = zeros(200,200);
y = zeros(200,200);
%���ƽ�壬����
for i = 1:200
    for j =1:100
        x(i,j) = ( xl(j) + xl(j+1) ) / 2;
        y(i,j) = ( yl(i) + yl(i+1) ) / 2;
    end
end
%�Ҳྲ��
for i = 1:200
    for j = 101:200
        x(i,j) = ( xr(j-100) + xr(j-100+1) ) / 2;
        y(i,j) = ( yr(i) + yr(i+1) ) / 2;
    end
end

%���ӵ��ٶ�
vx = zeros(200,200);
vy = zeros(200,200);
%����500m/s
for i = 1:200
    for j =1:100
        vx(i,j) = 500;
    end
end
%��ʼ����������
rho0 = 18.087e+3;
m = rho0 * 0.0001^2;
rho = zeros(200,200);
for i = 1:200
    for j = 1:200
        rho(i,j) = rho0;
    end
end

e = zeros(200,200);
for i = 1:200
    for j = 1:200
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
T = 201;
dt = 1e-7;
h = 0.0001*10.5;
for t = 1:T
    P1 = P;
    x1 = x;
    y1 = y;
    rho1 = rho;
    vx1 = vx;
    vy1 = vy;
    %����x�ı仯
    x = x + vx*dt;
    y = y + vy*dt;
    for iy = 1:200
        for ix = 1:200
            drho = 0;
            de = 0;
            dvx = 0;
            dvy = 0;
            for jy = 1:200
                for jx = 1:200
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
    %��e
   
    %��PH;c0 = 0.251;s = 1.51;
    PH = (c0^2*(1/rho0 - 1./rho))./(1/rho0 - s*(1/rho0 - 1./rho)).^2;
    %��eH
    eH = 0.5*PH.*(1/rho0 - 1./rho);
    %��P
    P = PH + Gamma.*rho.*(e - eH);
end
subplot(2,2,2);
plot(x,y,'o');
xlabel('x');
ylabel('y');
subplot(2,2,3);
plot(x(100,:),P(100,:));
xlabel('x');
ylabel('P(-sigma)');
subplot(2,2,4);
plot(x(100,:),vx(100,:));
xlabel('x');
ylabel('vx');