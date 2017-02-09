clear all; close all; clc;




Ecm_pair=linspace(1.05/0.511,3/0.511,100);
cos_theta_grid=linspace(-1,1,100);
qq=0;cross=[];X=[];Y=[];za=[];
% for cos_theta_cm=cos_theta_grid
for qq=1:1:length(cos_theta_grid)
%     qq=qq+1;
    cos_theta_cm=cos_theta_grid(qq);
   
    cross(:,qq)=breit_wheeler_crossection_in_the_CoM(Ecm_pair,cos_theta_cm);
    
        r0=2.818e-15;%[m]
%      % % \epsilon is the energy of each particle in CoM system in units of mc^2
%         epsilon=Ecm_pair./2;
%         % % \beta is understood to be the velocity of each colliding particle in the CoM system
%         beta=sqrt(1-1./epsilon.^2);
%         
%         Vq(:,qq)=((r0.^2.*beta)./(4.*epsilon.^2)).*...
%             ((1+2.*beta.^2.*(1-cos_theta_cm.^2)-beta.^4-beta.^4.*(1-cos_theta_cm.^2).^2)/(1-beta.^2.*cos_theta_cm.^2).^2);
%     
            
        Vq(:,qq)=((r0.^2.*sqrt(1-1./(Ecm_pair./2).^2))./(4.*(Ecm_pair./2).^2)).*...
            ((1+2.*(1-1./(Ecm_pair./2).^2).*(1-cos_theta_cm.^2)-(1-1./(Ecm_pair./2).^2).^2-(1-1./(Ecm_pair./2).^2).^2.*(1-cos_theta_cm.^2).^2)./(1-(1-1./(Ecm_pair./2).^2).*cos_theta_cm.^2).^2);
    


end

[X,Y]=meshgrid(Ecm_pair,cos_theta_grid);
figure(1)
mesh(X,Y,Vq')
% pcolor(Vq)
figure(2)
mesh(X,Y,cross')
      