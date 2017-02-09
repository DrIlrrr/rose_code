clear all; close all; clc;
%   stop
ifig=1;
for wp=4%:6;
    
    load(['data_WP_' num2str(wp) '.dat'],'-mat')
    
    
    lim_th=find(theta_max>0);
    
    theta_max_1=theta_max(lim_th);
    theta_max=[];
    theta_max=theta_max_1;
    bandwith_cm_1=bandwith_cm(lim_th);
    bandwith_cm=[];
    bandwith_cm=bandwith_cm_1;
    
  
    h=4.135667662e-15;% Plank constant [eV*s].
    c=3e8;% speed of light [m/s]
    hnu=h*c/laserwl;%laser energy in eV
    mc2=0.511e6;%electron rest mass in eV
    
    d1=4*gamma*hnu/mc2;
    
    T1=(((gamma.*theta_max).^2./sqrt(12))./(1+(gamma.*theta_max).^2./2)...
        +(norm_em_x^2/sigma_0^2)./(1+sqrt(12).*norm_em_x^2/sigma_0^2)...
        +(norm_em_y^2/sigma_0y^2)./(1+sqrt(12).*norm_em_y^2/sigma_0y^2)).*(1/(1+d1));
    
    T2=((2+d1)/(1+d1))*(delta_gamma/gamma);
    
    T3=(1/(1+d1))*laser_bandwidth*0;
    
    T4=((M2*laserwl)/(2*pi*W0))^2*0;
    
    T5=((a_0^2/3)/(1+a_0^2/2))*0;
    
    
    formBW=sqrt(T1.^2+T2^2+T3^2+T4^2+T5^2);
        
    figure(ifig)
    ifig=ifig+1;
    set(gca,'FontSize',16)
    hold on
    plot(theta_max,bandwith_cm,'-ok')
    %     plot(theta_max,New_my2,'--xb')
    plot(theta_max,formBW,'-sg')
    hold off
    grid on
    ylim([0 max([max(formBW) max(formBW)])])
    xlim([0 max(theta_max)])
    xlabel('\theta_{max}')
    ylabel('bandwidth')
    set(gca,'FontSize',14)
    legend('Simulation',...
        'Formula',...       
        'Location','NorthWest')
    title({['emitx=' num2str(e_emx_sqrt) '; emity=' num2str(e_emy_sqrt) ';' ],...
        ['Ee_{max}=' num2str(max(e_E)/1e6) ' [MeV]; \sigma_E=' num2str(std(e_E)/mean(e_E)) '; \nu_{max}=' num2str(max(full_spectrum)/1e3) ' [MeV];']})
    filename = ['NOENS_CE_' num2str(WP) ];
    fname = [ filename '.png'];
    print('-dpng', fname);
    
    
end
