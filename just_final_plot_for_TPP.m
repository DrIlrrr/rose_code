clear all; close all; clc;
rflags.PLOTS =1;
ifig=0;
% stop
save_dir=[pwd '/TPP_260_sigma_end_cut_1_L_0/nbin_x_21_y_21_z_21/'];
filename=[save_dir 'N_pairvents_'];
load([save_dir 'main_out_put.dat'],'-mat');
ad=load([save_dir 'out_put.dat'],'-mat');
% pair_info=[Vq E_3 E_4 theta_3 theta_4 phi_3 phi_4 cos_alpha cos_theta Ecm_pair gamma_cm_pair];




%% Make a histogram of E_CoM
n_bin=20;
ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
hist(ad.Ecm_tot,n_bin)
xlabel('E_{CoM}','FontSize',20)
ylabel('Number of pair','FontSize',20)
filename = [save_dir 'just_Final_fig_' num2str(ifig)];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
%     print('-depsc', fname);
print('-r300','-dpng', fname2);
%% and extract from it coordinates of center of bin

[N_pair_per_bin cc_bin]=hist(ad.Ecm_tot,n_bin);

N_pair=2*N_pair_per_bin*ad.weight_1*ad.weight_2;


%%



%     bin_windnes=1;
%     nbin_plot=floor((max(ad.Ecm_tot)*1e-6-min(ad.Ecm_tot)*1e-6)/bin_windnes);% nbin to make spectrum per KeV
%     xs_tot=linspace(min(ad.Ecm_tot)*1e-6,max(ad.Ecm_tot)*1e-6,nbin_plot);
%     %      bar(xs_tot,hist(ad.Ecm_tot,nbin_plot)*ad.weight_1*ad.weight_2,'grouped')%,'hist','g')
%
%     N_pair=hist(ad.Ecm_tot,nbin_plot)*ad.weight_1*ad.weight_2;
%
%     ifig=ifig+1;
%     if rflags.PLOTS ==1;
%         figure(ifig)
%     else
%         figure('visible','off');
%     end
%     hold on
%     bar(xs_tot,N_pair)%,'hist','g')
%     plot(xs_tot,N_pair,'-o')
%     hold off
%     set(gca,'FontSize',16)
%     grid on
%     xlim([0 max(xs_tot+10)])
%     xlabel('E_{CM}')
%     ylabel(['Number per ' num2str(bin_windnes) ' MeV'])
%     filename = [save_dir 'just_Final_fig_' num2str(ifig)];
%     fname = [ filename '.eps'];fname2 = [ filename '.png'];
%     %     print('-depsc', fname);
%     print('-r300','-dpng', fname2);



%     stop
xs_tot=cc_bin/1e6;
K=((xs_tot).^2-(0.511)^2)./(2.*(0.511)^2);




ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
set(gca,'FontSize',16)
hold on
%     bar(xs_tot,K)
plot(xs_tot,K,'-x')
hold off
grid on
xlim([0 max(xs_tot+10)])
xlabel('E_{CM}')
ylabel('K')
filename = [save_dir 'just_Final_fig_' num2str(ifig)];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
%     print('-depsc', fname);
print('-r300','-dpng', fname2);

%     if K>=4 && K<=4.6 % tut nado ispolsovat find ili cikl
%     f_1=(5.6+20.4.*(K-4)-10.9.*(K-4).^2-3.6.*(K-4).^3+7.4.*(K-4).^4).*1e-3.*(K-4).^2;
%     end
%     if (K > 4.6) && (K <= 6)
%     f_2=0.582814-0.29842.*K+0.04354.*K.^2-0.0012977.*K.^3;
%     end
%     if (K > 6) && (K <= 18)
%     f_3=(3.1247-1.3394.*K+0.14612.*K.^2)./(1+0.4648.*K+0.016683.*K.^2);
%     end
%     if K > 18
%     f_4=(28/9).*log(2.*K)-(218/27)+(1./K).*(-(4/3).*(log(2.*K)).^3+3.863.*(log(2.*K)).^2-11.*(log(2.*K))+27.9)
%     end
%
%
%
%     stop
%
%     sec_BH=14.*(3.11.*log(2.*K)-8.07);
%     sec_B=(14./K).*((4/3).*log(2.*K).^3-3.*log(2.*K).^2+6.84.*log(2.*K)-21.51);
%
%
%     % tottal_cross=1e-33.*(sec_BH-sec_B).*ad.weight_1.*ad.weight_2./(4.*pi.*delta_x.*delta_y.*1e4);
%     tottal_cross=1e-33.*(sec_BH-sec_B)./(4.*pi.*delta_x.*delta_y.*1e4);

%% cross-section block


alpha=1/137;
re=2.817e-15;

K=((xs_tot).^2-(0.511)^2)./(2.*(0.511)^2);
% K=[4:0.01:1000];
part_1=find(K>4 & K<=4.6);
f_1=(5.6+20.4.*(K(part_1)-4)-10.9.*(K(part_1)-4).^2-3.6.*(K(part_1)-4).^3+7.4.*(K(part_1)-4).^4).*1e-3.*(K(part_1)-4).^2;

part_2=find(K>4.6 & K<=6);
f_2=0.582814-0.29842.*K(part_2)+0.04354.*K(part_2).^2-0.0012977.*K(part_2).^3;

part_3=find(K>6 & K<=18);
f_3=(3.1247-1.3394.*K(part_3)+0.14612.*K(part_3).^2)./(1+0.4648.*K(part_3)+0.016683.*K(part_3).^2);

part_4=find(K>18);
f_4=(28/9).*log(2.*K(part_4))-(218/27)+(1./K(part_4)).*(-(4/3).*(log(2.*K(part_4))).^3+3.863.*(log(2.*K(part_4))).^2-11.*(log(2.*K(part_4)))+27.9);


tottal_cross=[re^2.*alpha.*f_1 re^2.*alpha.*f_2 re^2.*alpha.*f_3 re^2.*alpha.*f_4];
K_p=[K(part_1) K(part_2) K(part_3) K(part_4)];
X_p=[xs_tot(part_1) xs_tot(part_2) xs_tot(part_3) xs_tot(part_4)];

%     ifig=ifig+1;
%     if rflags.PLOTS ==1;
%         figure(ifig)
%     else
%         figure('visible','off');
%     end
%     ax.FontSize=20;
%     plot(K_p,tottal_cross,'LineWidth',2)
%     xlabel('K')
%     ylabel('F')
%     filename = [save_dir 'just_Final_fig_' num2str(ifig)];
%     fname = [ filename '.eps'];fname2 = [ filename '.png'];
%     %     print('-depsc', fname);
%     print('-r300','-dpng', fname2);
%
%
%     ifig=ifig+1;
%     if rflags.PLOTS ==1;
%         figure(ifig)
%     else
%         figure('visible','off');
%     end
%     plot(X_p,tottal_cross,'LineWidth',2)
%     xlabel('E_{Com}')
%     ylabel('F')
%     filename = [save_dir 'just_Final_fig_' num2str(ifig)];
%     fname = [ filename '.eps'];fname2 = [ filename '.png'];
%     %     print('-depsc', fname);
%     print('-r300','-dpng', fname2);
%
%
%     ifig=ifig+1;
%     if rflags.PLOTS ==1;
%         figure(ifig)
%     else
%         figure('visible','off');
%     end
%     plot(X_p,K_p,'LineWidth',2)
%     xlabel('E_{Com}')
%     ylabel('K')
%     filename = [save_dir 'just_Final_fig_' num2str(ifig)];
%     fname = [ filename '.eps'];fname2 = [ filename '.png'];
%     %     print('-depsc', fname);
%     print('-r300','-dpng', fname2);
%%
ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
set(gca,'FontSize',16)
hold on
%     plot(K_p,tottal_cross,'LineWidth',2)
bar(K_p,tottal_cross)
plot(K_p,tottal_cross,'-x')
hold off
grid on
%     xlim([0 max(xs_tot+10)])
xlabel('K')
ylabel('Tot cross section')
filename = [save_dir 'just_Final_fig_' num2str(ifig)];
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
bar(xs_tot,tottal_cross)
plot(xs_tot,tottal_cross,'-x')
hold off
grid on
xlim([0 max(xs_tot+10)])
xlabel('E_{CM}')
ylabel('Tot cross section')
filename = [save_dir 'just_Final_fig_' num2str(ifig)];
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
bar(xs_tot,tottal_cross.*N_pair./(ad.delta_x.*ad.delta_y))
plot(xs_tot,tottal_cross.*N_pair./(ad.delta_x.*ad.delta_y),'-o')
hold off
grid on
xlim([0 max(xs_tot+10)])
title(['sum ' num2str(sum(tottal_cross.*N_pair./(ad.delta_x.*ad.delta_y))) ])
xlabel('E_{CM}')
ylabel('Tot cross * num')
filename = [save_dir 'just_Final_fig_' num2str(ifig)];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
%     print('-depsc', fname);
print('-r300','-dpng', fname2);


%%
ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end

nbin_plot=20;
xs_gtot=linspace(min(ad.gammacm_tot),max(ad.gammacm_tot),nbin_plot);
bar(xs_gtot,hist(ad.gammacm_tot,nbin_plot)*ad.weight_1*ad.weight_2,'grouped')%,'hist','g')
grid on
set(gca,'FontSize',16)
xlabel('\gamma_{cm} total')
xlim([0 max(xs_gtot)+0.1*max(xs_gtot)])
% ylabel(['Number per ' num2str(bin_windnes) ' MeV'])
filename = [save_dir 'just_Final_fig_' num2str(ifig)];
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
xs_gtot=linspace(min(ad.gammacm_tot),max(ad.gammacm_tot),nbin_plot);
bar(xs_gtot,hist(ad.gammacm_tot,nbin_plot)*ad.weight_1*ad.weight_2,'grouped')%,'hist','g')
grid on
set(gca,'FontSize',16)
xlabel('\gamma_{cm} total')
xlim([0 max(xs_gtot)+0.1*max(xs_gtot)])
ylabel(['N_{pair}'])
subplot 312
hold on
bar(xs_tot,N_pair)%,'hist','g')
plot(xs_tot,N_pair,'-o')
hold off
set(gca,'FontSize',16)
grid on
xlim([0 max(xs_tot+10)])
xlabel('E_{CM}')
ylabel(['N_{pair}'])

subplot 313
set(gca,'FontSize',16)
hold on
bar(xs_tot,tottal_cross.*N_pair)
plot(xs_tot,tottal_cross.*N_pair,'-o')
hold off
grid on
xlim([0 max(xs_tot+10)])
title(['sum ' num2str(sum(tottal_cross.*N_pair./(ad.delta_x.*ad.delta_y))) ])
xlabel('E_{CM}')
ylabel('N event')
filename = [save_dir 'just_Final_fig_' num2str(ifig)];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
%     print('-depsc', fname);
print('-r300','-dpng', fname2);
