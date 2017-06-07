clear all; close all; clc;
global ev cs cross


 ev=load('e_grid.dat','ascii')';
 cs=load('w_grid.dat','ascii')';
 cross=load('cs_grid.dat','ascii')';

%% 
ev_m=linspace(0.3,2,60);
cthv_m=linspace(-0.9995,0.9995,60);
 

[Eq,Wq] = meshgrid(ev_m,cthv_m);


% Eq=ev;
% Wq=cs;

Vq=interp2(ev,cs,cross,Eq,Wq,'linear');



figure(90)
subplot(2,2,1)
surf(Wq,Eq,Vq)
title('interpol')
subplot(2,2,2)
 mesh(cs(:,1),ev(1,:),cross')
 title('data')
subplot(2,2,[3 4])
hold on
surf(Wq,Eq,Vq)
%  scatter3(Wq,Eq,Vq,'r')
% scatter3(grid300(:,1),grid300(:,2),grid300(:,3),'.b')
 mesh(cs(:,1),ev(1,:),cross')%,'EdgeColor','none')
hold off
xlabel('Cos(\theta)')
ylabel('E_{cm}')
view(44,30)

% stop
E_cm=1.9;
for ni=1:1:5000
[theta_out]=rejection_method(E_cm);

th(ni)=theta_out;
end


figure(100)
hist(th,20)


% 
% 
% interpn(grid300(:,1),acos(grid300(:,2)),grid300(:,3),Xq,Yq)



