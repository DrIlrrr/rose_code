% function plot_comton_cross_com
clear all; close all; clc;
ifig=1;


% E_grid=linspace(0,2.4,100);
% 
% % E_grid=linspace(2.0001,1e2,100);
% cos_theta_grid=linspace(-1,1,100);


ev=load('rs_grid.dat','ascii')';
cs=load('w_grid.dat','ascii')';
crosst=load('cs_grid.dat','ascii')';

E_grid=ev(1,:);
cos_theta_grid=cs(:,1);

qq=0;cross=[];X=[];Y=[];
for Ecm_pair=E_grid
    qq=qq+1;
    
    cross(qq,:)=aprox_gamma_gamma_crossection_in_the_CoM(Ecm_pair,cos_theta_grid);
    
end



ifig=ifig+1;
figure(ifig)
mesh((cs(:,1)),ev(1,:),1e-34.*crosst')
xlabel('\theta')
ylabel('E_{cm}')
view(44,30)
fname2 = [ 'cross_' num2str(ifig) '.png'];
print('-r300','-dpng', fname2);

[X,Y]=meshgrid(E_grid,cos_theta_grid);


ifig=ifig+1;
figure(ifig)
mesh(X,Y,cross')
ylabel('\theta')
xlabel('E_{cm}')

fname2 = [ 'cross_' num2str(ifig) '.png'];
print('-r300','-dpng', fname2);


[m,n]=size(cross);
tot_e_cross=[];
for ni=1:1:m
    tot_e_cross(ni)=(2*pi)*sum(cross(ni,:)).*(cos_theta_grid(2)-cos_theta_grid(1));
    tot_e_crossT(ni)=(2*pi)*sum(1e-34.*crosst(:,ni)).*(cos_theta_grid(2)-cos_theta_grid(1));
end
en_x=E_grid;


ifig=ifig+1;
figure(ifig)
% subplot 211
hold on
plot(E_grid,tot_e_cross,'-r','LineWidth',2)
plot(E_grid,tot_e_crossT,'--b','LineWidth',2)
hold off
legend('approx','table',0)
subplot 212
% plot(en_x,tot_e_cross-tot_e_crossT,'LineWidth',2)
plot(en_x,tot_e_cross,'-r','LineWidth',2)
plot(linspace(0.1,2.4,100),tot_e_crossT,'--b','LineWidth',2)
ylabel('cross-section','FontSize',20)
xlabel('E_{CoM}','FontSize',20)
fname2 = [ 'cross_' num2str(ifig) '.png'];
print('-r300','-dpng', fname2);



ifig=ifig+1;
figure(ifig)
% hold on
% plot(en_x,tot_e_cross,'LineWidth',2)
loglog(E_grid,tot_e_crossT,'LineWidth',2)
% hold off
% legend('approx','table',0)
ylabel('cross-section','FontSize',20)
xlabel('E_{CoM}','FontSize',20)
