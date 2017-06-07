function [cross]=aprox_gamma_gamma_crossection_in_the_CoM(Ecm_pair,cos_theta_cm)
% approximation of the differential cross-section 
%for unpolarized photon in the CoM


r0=2.818e-15;%[m]
alpha=1/137;
k=Ecm_pair./2;
cross=(alpha.^2.*r0.^2./(4.*pi.^2)).*(139/90^2).*k.^6.*(3+cos_theta_cm.^2).^2.*(1+(160/139).*(k.^2.*(1-cos_theta_cm.^2))./(3+cos_theta_cm));
% cross=(alpha.^2.*r0.^2./(4.*pi)).*(139/90^2).*k.^6.*(3+cos_theta_cm.^2).^2.*(1+(160/139).*(k.^2.*(1-cos_theta_cm.^2))./(3+cos_theta_cm));