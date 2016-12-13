clear all; close all; clc;

Ecm_pair_s=linspace(1.05/0.511,10/0.511,100);
cos_theta_cm_s=linspace(-1,1,100);


qn=0;

for ni=1:1:100
    qn=qn+1;
    qm=0;
    for mi=1:1:100
        qm=qm+1;
        Ecm_pair=Ecm_pair_s(qn);
        cos_theta_cm=cos_theta_cm_s(qm);
        
        
        %% version from formula tesi Davide chapter 3.2
%         sin_theta_cm_2=1-cos_theta_cm.^2;
        %  mc2=0.511e6;% electron rest mass [eV]
        r0=2.818e-15;%[m]
        
        
        % % \epsilon is the energy of each particle in CoM system in units of mc^2
        epsilon=Ecm_pair./2;
        % % \beta is understood to be the velocity of each colliding particle in the CoM system
        beta=sqrt(1-1./epsilon.^2);
        
        
        %         Vq(qn,qm)=((r0.^2.*beta)./(4.*epsilon.^2)).*...
        %             ((1+2.*beta.^2.*sin_theta_cm_2-beta.^4-beta.^4.*sin_theta_cm_2.^2)/(1-beta.^2.*cos_theta_cm.^2).^2);
        
        Vq(qn,qm)=((r0.^2.*beta)./(4.*epsilon.^2)).*...
            ((1+2.*beta.^2.*(1-cos_theta_cm.^2)-beta.^4-beta.^4.*(1-cos_theta_cm.^2).^2)./(1-beta.^2.*cos_theta_cm.^2).^2);
        
        cross(qn,qm)=breit_wheeler_crossection_in_the_CoM(Ecm_pair,cos_theta_cm);
        %za(qn,qm)= ((1+2.*(1-1./(Ecm_pair./2).^2).*(1-cos_theta_cm.^2)-(1-1./(Ecm_pair./2).^2).^2-(1-1./(Ecm_pair./2).^2).^2.*(1-cos_theta_cm.^2).^2)/(1-(1-1./(Ecm_pair./2).^2).*cos_theta_cm.^2).^2);
    end
end

[X,Y]=meshgrid(cos_theta_cm_s,Ecm_pair_s);
figure(1)
mesh(X,Y,Vq)
% pcolor(Vq)
figure(2)
mesh(X,Y,cross)
% figure(3)
% mesh(X,Y,za)

