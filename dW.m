function [dWx,dWy] = dW(Rx,Ry,h)% Rx = xi-xj Ry = yi-yj
%2D�㷨
%ÿ��С����Ϊdx,h=dx*n
ad = 5/(pi*h^2);%2D����
R = sqrt(Rx^2+Ry^2)/h;
% if R<=1 && R>0
    dWx = ad * (1-R)^2 * (-12*R) * Rx/(h*sqrt(Rx^2+Ry^2));
    dWy = ad * (1-R)^2 * (-12*R) * Ry/(h*sqrt(Rx^2+Ry^2));
% else
%     dWx = 0;
%     dWy = 0;
% end
end