function [cross]=compton_crossection_in_the_CoM(Ecm_pair,cos_theta_cm)

%% Compton crossection in the CoM
%  mc2=0.511e6;% electron rest mass [eV]
r0=2.818e-15;%[m]

%% version from formula
% K1=(-(Ecm_pair.^2-mc2.^2).*cos_theta_cm+(Ecm_pair.^4-mc2.^4))./(mc2.^2.*2.*Ecm_pair.^2);
%
% X=mc2^2./(Ecm_pair.^2-mc2^2);% X=1/K=(Ecm_pair^2-mc2^2)/mc2^2;
%
% Y=1./K1;
%
% X_grand=0.25.*(X./Y+Y./X)+(X-Y)+(X-Y).^2;
%
% cross=4.*pi.*(mc2)^2.*r0^2.*X_grand./Ecm_pair.^2;

%% version from fortran
mc2=(0.511e6)^2;% electron rest mass [eV]
s=Ecm_pair^2;
ekappa=-(mc2-s)./mc2;
ekprimo=(-(s-mc2).^2.*cos_theta_cm+(s.^2-mc2^2))./(mc2.*2.*s);
x=1./ekappa;
y=1./ekprimo;
x1=ekprimo./ekappa;
xgrande=(0.25*(x1+1./x1)+(x-y)+(x-y).^2);
cross=2.*mc2.*r0^2.*xgrande./s;

