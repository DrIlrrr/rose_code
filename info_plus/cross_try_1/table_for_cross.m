clear all; close all; clc;

 load('grid300.dat','ascii')
% scatter3(grid300(:,1),grid300(:,2),grid300(:,3))


cross=load('grid_linear.dat','ascii')';

%% 

V=cross;%magic(10);

ev=linspace(0.3,2,301);
cthv=linspace(-0.9995,0.9995,301)


[X,Y] = meshgrid(ev,cthv);


for ni=50000;%1:100;
%     cross(3,1)
Xq=grid300(ni,1)%linspace(0.3,2,10)
Yq=grid300(ni,2)%linspace(-0.9995,0.9995,10)
% zz(ni)=
grid300(ni,3)
Vq=interp2(X,Y,V,Xq,Yq,'linear')
% za(ni)=Vq
end



figure(90)
hold on
scatter3(Xq,Yq,Vq,'r')
scatter3(grid300(:,1),grid300(:,2),grid300(:,3),'.b')
surf(ev,cthv,V,'EdgeColor','none')
hold off
xlabel('x')
ylabel('y')
view(-35,10)



% 
% 
% interpn(grid300(:,1),acos(grid300(:,2)),grid300(:,3),Xq,Yq)



