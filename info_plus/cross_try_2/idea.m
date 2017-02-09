clear all; close all; clc;

V=magic(5)
%% 
x=linspace(-1,1,5);
y=linspace(-1,1,5);
%  cthv_m=linspace(-0.9995,0.9995,6)
 

 [X,Y] = meshgrid(x,y);


Xq=x(2)+0.28
Yq=y(2)+0.08
Vq=interp2(X,Y,V,Xq,Yq,'linear')


figure(90)
hold on
scatter3(Xq,Yq,Vq,'r')
surf(x,y,V,'EdgeColor','none')
hold off
xlabel('x')
ylabel('y')
