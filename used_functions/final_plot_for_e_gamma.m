function final_plot_for_e_gamma(pcm_tot,gammacm_tot,Ecm_tot,weight_1,weight_2,delta_x,delta_y)
% plot final plot for E_cm Gamma_cm_and N events
global ifig save_dir rflags
if rflags.e_gamma==1
    
    ifig=ifig+1;
    if rflags.PLOTS ==1;
        figure(ifig)
    else
        figure('visible','off');
    end
    
    nbin_plot=20;
    xs_gtot=linspace(min(gammacm_tot),max(gammacm_tot),nbin_plot);
    bar(xs_gtot,hist(gammacm_tot,nbin_plot)*weight_1*weight_2,'grouped')%,'hist','g')
    grid on
    set(gca,'FontSize',16)
    xlabel('\gamma_{cm} total')
    xlim([0 max(xs_gtot)+0.1*max(xs_gtot)])
    % ylabel(['Number per ' num2str(bin_windnes) ' MeV'])
    filename = [save_dir 'Final_fig_' num2str(ifig)];
    fname = [ filename '.eps'];fname2 = [ filename '.png'];
    %     print('-depsc', fname);
    print('-r300','-dpng', fname2);
    
    
    %% momentum plot
    
    % mom_cm=sqrt(gammacm_tot.^2-1);
    
    %
    % ifig=ifig+1;
    % if rflags.PLOTS ==1;
    %     figure(ifig)
    % else
    %     figure('visible','off');
    % end
    %
    %     nbin_plot=20;
    %     xs_gtot=linspace(min(pcm_tot),max(pcm_tot),nbin_plot);
    %     bar(xs_gtot,hist(pcm_tot,nbin_plot)*weight_1*weight_2,'grouped')%,'hist','g')
    %     grid on
    %     set(gca,'FontSize',16)
    %     xlabel('P_{cm} total')
    %     % xlim([0 max(xs_gtot)+0.1*max(xs_gtot)])
    %     % ylabel(['Number per ' num2str(bin_windnes) ' MeV'])
    %     filename = [save_dir 'Final_fig_' num2str(ifig)];
    %     fname = [ filename '.eps'];fname2 = [ filename '.png'];
    %     %     print('-depsc', fname);
    %     print('-r300','-dpng', fname2);
    
    %%
    
    
    
    bin_windnes=25;
    nbin_plot=floor((max(Ecm_tot)*1e-6-min(Ecm_tot)*1e-6)/bin_windnes);% nbin to make spectrum per KeV
    xs_tot=linspace(min(Ecm_tot)*1e-6,max(Ecm_tot)*1e-6,nbin_plot);
    % bar(xs_tot,hist(Ecm_tot,nbin_plot)*weight_1*weight_2,'grouped')%,'hist','g')
    
    
    num_e=hist(Ecm_tot,nbin_plot)*weight_1*weight_2;
    
    ifig=ifig+1;
    if rflags.PLOTS ==1;
        figure(ifig)
    else
        figure('visible','off');
    end
    hold on
    bar(xs_tot,num_e)%,'hist','g')
    plot(xs_tot,num_e,'-o')
    hold off
    set(gca,'FontSize',16)
    grid on
    xlim([0 max(xs_tot+10)])
    xlabel('E_{CM}')
    ylabel(['Number per ' num2str(bin_windnes) ' MeV'])
    filename = [save_dir 'Final_fig_' num2str(ifig)];
    fname = [ filename '.eps'];fname2 = [ filename '.png'];
    %     print('-depsc', fname);
    print('-r300','-dpng', fname2);
    % stop
    K=((xs_tot).^2)./(2.*0.511.*105.6);
    
    ifig=ifig+1;
    if rflags.PLOTS ==1;
        figure(ifig)
    else
        figure('visible','off');
    end
    set(gca,'FontSize',16)
    hold on
    bar(xs_tot,K)
    plot(xs_tot,K,'-x')
    hold off
    grid on
    xlim([0 max(xs_tot+10)])
    xlabel('E_{CM}')
    ylabel('K')
    filename = [save_dir 'Final_fig_' num2str(ifig)];
    fname = [ filename '.eps'];fname2 = [ filename '.png'];
    %     print('-depsc', fname);
    print('-r300','-dpng', fname2);
    
    
    sec_BH=14.*(3.11.*log(2.*K)-8.07);
    sec_B=(14./K).*((4/3).*log(2.*K).^3-3.*log(2.*K).^2+6.84.*log(2.*K)-21.51);
    
    
    % tottal_cross=1e-33.*(sec_BH-sec_B).*weight_1.*weight_2./(4.*pi.*delta_x.*delta_y.*1e4);
    tottal_cross=1e-33.*(sec_BH-sec_B)./(4.*pi.*delta_x.*delta_y.*1e4);
    
    ifig=ifig+1;
    if rflags.PLOTS ==1;
        figure(ifig)
    else
        figure('visible','off');
    end
    set(gca,'FontSize',16)
    hold on
    bar(xs_tot,tottal_cross)
    plot(xs_tot,tottal_cross,'-x')
    hold off
    grid on
    xlim([0 max(xs_tot+10)])
    xlabel('E_{CM}')
    ylabel('Tot cross section')
    filename = [save_dir 'Final_fig_' num2str(ifig)];
    fname = [ filename '.eps'];fname2 = [ filename '.png'];
    %     print('-depsc', fname);
    print('-r300','-dpng', fname2);
    
    
    ifig=ifig+1;
    if rflags.PLOTS ==1;
        figure(ifig)
    else
        figure('visible','off');
    end
    set(gca,'FontSize',16)
    hold on
    bar(xs_tot,tottal_cross.*num_e)
    plot(xs_tot,tottal_cross.*num_e,'-o')
    hold off
    grid on
    xlim([0 max(xs_tot+10)])
    title(['sum ' num2str(sum(tottal_cross.*num_e)) ])
    xlabel('E_{CM}')
    ylabel('Tot cross * num')
    filename = [save_dir 'Final_fig_' num2str(ifig)];
    fname = [ filename '.eps'];fname2 = [ filename '.png'];
    %     print('-depsc', fname);
    print('-r300','-dpng', fname2);
    
    
    
    
    ifig=ifig+1;
    if rflags.PLOTS ==1;
        figure(ifig)
    else
        figure('visible','off');
    end
    subplot 311
    
    nbin_plot=20;
    xs_gtot=linspace(min(gammacm_tot),max(gammacm_tot),nbin_plot);
    bar(xs_gtot,hist(gammacm_tot,nbin_plot)*weight_1*weight_2,'grouped')%,'hist','g')
    grid on
    set(gca,'FontSize',16)
    xlabel('\gamma_{cm} total')
    xlim([0 max(xs_gtot)+0.1*max(xs_gtot)])
    ylabel(['N_{pair}'])
    subplot 312
    hold on
    bar(xs_tot,num_e)%,'hist','g')
    plot(xs_tot,num_e,'-o')
    hold off
    set(gca,'FontSize',16)
    grid on
    xlim([0 max(xs_tot+10)])
    xlabel('E_{CM}')
    ylabel(['N_{pair}'])
    
    subplot 313
    set(gca,'FontSize',16)
    hold on
    bar(xs_tot,tottal_cross.*num_e)
    plot(xs_tot,tottal_cross.*num_e,'-o')
    hold off
    grid on
    xlim([0 max(xs_tot+10)])
    title(['sum ' num2str(sum(tottal_cross.*num_e)) ])
    xlabel('E_{CM}')
    ylabel('N event')
    filename = [save_dir 'Final_fig_' num2str(ifig)];
    fname = [ filename '.eps'];fname2 = [ filename '.png'];
    %     print('-depsc', fname);
    print('-r300','-dpng', fname2);
    
    
end

