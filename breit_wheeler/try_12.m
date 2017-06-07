
clear all; close all; clc;
Ecm_pair=linspace(1.05/0.511,3/0.511,10);

cos_theta_cm=linspace(-1,1,10);


%  mc2=0.511e6;% electron rest mass [eV]
r0=2.818e-15;%[m]

%% version from formula tesi Davide chapter 3.2
sin_theta_cm_2=1-cos_theta_cm.^2;

% \epsilon is the energy of each particle in CoM system in units of mc^2
epsilon=Ecm_pair./2;
% \beta is understood to be the velocity of each colliding particle in the CoM system
beta=sqrt(1-1./epsilon.^2);

part_1=((r0.^2.*beta)./(4.*epsilon.^2));
Vq(:,:)=part_1.*...
    ((1+2.*beta.^2.*sin_theta_cm_2-beta.^4-beta.^4.*sin_theta_cm_2.^2)./(1-beta.^2.*cos_theta_cm.^2).^2);




















% ((r0.^2.*beta)./(4.*epsilon.^2)).*...
%             ((1+2.*beta.^2.*(1-cos_theta_cm.^2)-beta.^4-beta.^4.*(1-cos_theta_cm.^2).^2)/(1-beta.^2.*cos_theta_cm.^2).^2)






