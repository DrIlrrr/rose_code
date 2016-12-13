function [cross]=moller_crossection_in_the_CoM(Ecm_pair,cos_theta_cm)
global gvar
% cos_theta_cm
% stop
% if cos_theta_cm>0
%     stop
%     r0=0;%[m]
% else
    r0=2.818e-15;%[m]


%% version from formula
mc2=(0.511e6);% electron rest mass [eV]
%s=Ecm_pair^2;
% p=4*Ecm_pair
gamma_c=Ecm_pair./(2*mc2);
betta=sqrt(1-1./gamma_c.^2);
cross=r0^2.*((1+betta.^2)./(4.*betta.^4.*gamma_c.^2)).*...
    (4./(1-cos_theta_cm.^2).^2-3./(1-cos_theta_cm.^2)+(betta.^2./(1+betta.^2)).*(1+4./sqrt(1-cos_theta_cm.^2)));




  indc1=find(abs(cos_theta_cm)>gvar.lim_for_cos_theta_cm);
  cross(:,indc1)=0;
%  indc2=find(cos_theta_cm>0);
%  cross(:,indc2)=0;