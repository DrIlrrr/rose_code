% function plot_comton_cross_com
clear all; close all; clc;
E_grid=linspace(1.022e6,6e6,100);

% E_grid=linspace(2.0001,1e2,100);
cos_theta_grid=linspace(-1,1,100);
qq=0;cross=[];X=[];Y=[];
for Ecm_pair=E_grid
    qq=qq+1;
    
   
    cross(qq,:)=breit_wheeler_crossection_in_the_CoM(Ecm_pair,cos_theta_grid);
    
end

[X,Y]=meshgrid(E_grid,cos_theta_grid);
figure(1)
mesh(X,Y,cross')




[m,n]=size(cross);
tot_e_cross=[];
for ni=1:1:n
    tot_e_cross(ni)=(2*pi)*sum(cross(ni,:)).*(cos_theta_grid(2)-cos_theta_grid(1));
end
en_x=E_grid;


figure(2)

plot(en_x,tot_e_cross,'LineWidth',2)
ylabel('cross-section','FontSize',20)
xlabel('E_{CoM}','FontSize',20)







% 
% 
% 
% [m,n]=size(cross);
% tot_e_cross=[];
% for ni=1:1:n
%     tot_e_cross(ni)=(2*pi)*sum(cross(ni,:)).*(cos_theta_grid(2)-cos_theta_grid(1));
% end
% en_x=E_grid;
% 
% 
% figure(2)
% 
% plot(en_x,tot_e_cross,'LineWidth',2)
% ylabel('cross-section','FontSize',20)
% xlabel('E_{CoM}','FontSize',20)
% 
% 
% 
% 
% 
% % mesh(cross)
% 
% %
% %
% 
% %
% %  plot(Y)