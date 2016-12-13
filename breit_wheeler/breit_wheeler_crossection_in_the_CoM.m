function [Vq]=breit_wheeler_crossection_in_the_CoM(Ecm_pair,cos_theta_cm)

%% Compton crossection in the CoM
%  mc2=0.511e6;% electron rest mass [eV]
r0=2.818e-15;%[m]

% %% version from formula tesi Davide chapter 3.2
% sin_theta_cm_2=1-cos_theta_cm.^2;
% 
% % \epsilon is the energy of each particle in CoM system in units of mc^2
% epsilon=Ecm_pair./2;
% % \beta is understood to be the velocity of each colliding particle in the CoM system
% beta=sqrt(1-1./epsilon.^2);
% 
% 
% Vq=((r0.^2.*beta)./(4.*epsilon.^2)).*...
%     ((1+2.*beta.^2.*sin_theta_cm_2-beta.^4-beta.^4.*sin_theta_cm_2.^2)/(1-beta.^2.*cos_theta_cm.^2).^2);


%%
% \epsilon is the energy of each particle in CoM system in units of mc^2
epsilon=(Ecm_pair/0.511e6)./2;
% \beta is understood to be the velocity of each colliding particle in the CoM system
betta=sqrt(1-1./epsilon.^2);%beta is function in matlab use betta


Vq=((r0.^2.*betta)./(4.*epsilon.^2)).*...
    ((1+2.*betta.^2.*(1-cos_theta_cm.^2)-betta.^4-betta.^4.*(1-cos_theta_cm.^2).^2)./(1-betta.^2.*cos_theta_cm.^2).^2);