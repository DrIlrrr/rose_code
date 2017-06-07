clear all; close all; clc;
Ecm_pair=linspace(1.05/0.511,3/0.511,100);

% E_grid=linspace(2.0001,1e2,100);
cos_theta_cm=linspace(-1,1,100);
%% Compton crossection in the CoM
%  mc2=0.511e6;% electron rest mass [eV]
r0=2.818e-15;%[m]

%% version from formula tesi Davide chapter 3.2
sin_theta_cm_2=1-cos_theta_cm.^2;

% \epsilon is the energy of each particle in CoM system in units of mc^2
epsilon=1;%Ecm_pair./2;
% \beta is understood to be the velocity of each colliding particle in the CoM system
% beta=sqrt(1-1./epsilon.^2);

qq=0;
 beta1=linspace(0,1,100);
for beta=beta1
    qq=qq+1;
    
Vq(qq,:)=((r0.^2.*beta)./(4.*epsilon.^2)).*...
    ((1+2.*beta.^2.*sin_theta_cm_2-beta.^4-beta.^4.*sin_theta_cm_2.^2)/(1-beta.^2.*cos_theta_cm.^2).^2);
end
% stop

[X,Y]=meshgrid(beta1,cos_theta_cm);
figure(1)
mesh(X,Y,Vq)