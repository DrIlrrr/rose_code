% one beam several events
clear all; close all; clc;
make_path
%   stop
%for new global
global home_dir DIRECTORY_FOR_CAIN BASE_DIRECTORY;
global rflags
[rflags] = flags_for_run;
home_dir=[pwd '/CAIN/'];
rflags.PLOTS =0;

ifig=1;
rflags.PLOTS =1;
just_plots=1;
qq=0;
scan_spectrum=[];
phot_angle=[];


aaa=[0:0.5:5];
rflags.pulseE=0.2;
scan_var_name=['\alpha_0'];
for var_for_scan=aaa
    qq=qq+1;
    rflags.angle=var_for_scan*(pi/180);
    %rflags.laserwl=515; % laser wavelenth [nm] nano meters
%    rflags.chargebunch = 250e-12;%Charge per electrons bunch [c]
%    defoc_param=1;%defocusing parameter by defoult 1 is no defocusing
     BASE_DIRECTORY = [pwd '/Beam_l_Echarge_' num2str(rflags.chargebunch*1e12) '_angle_' num2str(rflags.angle) '/'];
    
    
    [full_spectrum,phot_angle,weigth,x_phot,y_phot,z_phot,xp_phot,yp_phot,zp_phot,Sx,Sy,Sz]=scan_photons_plots(BASE_DIRECTORY);
    
    eval(['full_spectrum_' int2str(qq) '=full_spectrum;'])
    eval(['phot_angle_' int2str(qq) '=phot_angle;'])
    eval(['x_coor' int2str(qq) '=x_phot;'])
    eval(['y_coor' int2str(qq) '=y_phot;'])
    eval(['z_coor' int2str(qq) '=z_phot;'])
    eval(['xp_' int2str(qq) '=xp_phot;'])
    eval(['yp_' int2str(qq) '=yp_phot;'])
    eval(['zp_' int2str(qq) '=zp_phot;'])
    eval(['Sx_' int2str(qq) '=Sx;'])
    eval(['Sy_' int2str(qq) '=Sy;'])
    eval(['Sz_' int2str(qq) '=Sz;'])
    eval(['weigth_' int2str(qq) '=weigth;'])
    
%% Theor flux
    %  1  2         3     4    5    6    7     8      9        10       11    12 13 14
    %  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss
    electron_data=load([BASE_DIRECTORY 'cain_tmp/exp.dat']);
    std_x_size=std(electron_data(:,5));
    std_y_size=std(electron_data(:,6));
    std_s=std(electron_data(:,4));
    mean_x_el(qq)=mean(electron_data(:,5));
    mean_y_el(qq)=mean(electron_data(:,6));
    mean_s_el(qq)=mean(electron_data(:,4));
    number_electrons=electron_data(1,3)*length(electron_data(:,3));

    
    SPEED_OF_LIGHT=3e8;
    Sigma_th=0.665e-28;%[m^2] Thomsom cross section
    h=2*pi*1.054e-34;%Planc const [J*s]
    lambda_l=1029e-9;% [m] laser wave length
    pulseE=0.2; %[J]
    photons_number=pulseE/((h*SPEED_OF_LIGHT)/lambda_l);
    
    angle=var_for_scan;
    sigLr=14e-6;
    %         laser_length=var_val;
    sigt=2.5e-12;
    laser_length=sigt*SPEED_OF_LIGHT;
    
    %     delta_s=0;
    %     delta_x(qq)=rflags.shifting_laser_x;
    %     %     Ax(qq)=exp((delta_x(qq).^2)/(-2*(std_y_size.^2+sigLr^2)));
    As(qq)=exp(((mean_s_el(qq)).^2.*tan((1/2)*angle.*(pi/180))^2)/(-2*(std_x_size.^2+sigLr^2+((std_s)^2+(sigt*SPEED_OF_LIGHT)^2)*tan((1/2)*angle*(pi/180))^2)));    
    Ay(qq)=exp(((mean_y_el(qq)).^2)/(-2*(std_x_size.^2+sigLr^2)));
    Ax(qq)=exp(((mean_x_el(qq)).^2)/(-2*(std_x_size.^2+sigLr^2+((std_s)^2+(sigt*SPEED_OF_LIGHT)^2)*tan((1/2)*angle*(pi/180))^2)));
    NUM_ph_x(qq)=As(qq).*Ax(qq).*Ay(qq).*((Sigma_th*photons_number*number_electrons)/(2*pi*sqrt(std_y_size^2+sigLr^2)))...
        /sqrt((std_x_size^2+sigLr^2)+((std_s)^2+(sigt*SPEED_OF_LIGHT)^2)*tan((1/2)*angle*(pi/180))^2);
    %%
    scan_var(qq)=var_for_scan;%std_x_size(qq);
end



%% plot total pectrum
figure(ifig)
ifig=ifig+1;
cc1=jet(qq);
for ni=1:1:qq
    
    hold on
    plot(linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))),...
        histc(eval(['full_spectrum_' int2str(ni) ]),linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))))*eval(['weigth_' int2str(ni) ]),...
        'color',cc1(ni,:),'LineWidth',2)
    
    % plot(linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))),...
    %     histc(eval(['full_spectrum_' int2str(ni) ]),linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))))*weigth,...
    %     'LineWidth',2);
    
    grid on
    hold off
end
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')
legendCell = cellstr(num2str(aaa', 'def=%-f'))
legend(legendCell,'FontSize',16)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);
ni=[];

%%
ww=0;
for ni=aaa
    ww=ww+1
    eval(['aa_' num2str(ww) '=find(abs(phot_angle_' num2str(ww) '));']);
    total_num_in_b(ww)=eval(['length(full_spectrum_' int2str(ww) '(aa_' int2str(ww) '))*weigth_' int2str(ww) ';']);
    
end

theta=0.0043;
ww=0;
for ni=aaa
    ww=ww+1
    eval(['aa_' num2str(ww) '=find(abs(phot_angle_' num2str(ww) ')<' num2str(theta) ');']);
    num_in_b(ww)=eval(['length(full_spectrum_' int2str(ww) '(aa_' int2str(ww) '))*weigth_' int2str(ww) ';']);
    
end


figure(ifig)
ifig=ifig+1;
subplot 211
hold on
plot(scan_var,total_num_in_b,'-.xr','LineWidth',0.5)
plot(scan_var,total_num_in_b,'xr','LineWidth',3)
hold off
grid on
title('total')
ylim([0 max(total_num_in_b)+max(total_num_in_b)*1e-1])
set(gca,'FontSize',16)
ylabel('Flux')
% ylabel('number of scattered photons')
subplot 212
hold on
plot(scan_var,num_in_b,'-.xb','LineWidth',0.5)
plot(scan_var,num_in_b,'xb','LineWidth',3)
title(['in Theta=' num2str(theta) ])
hold off
grid on
ylim([0 max(num_in_b)+max(num_in_b)*1e-1])
set(gca,'FontSize',16)
ylabel('Flux')
% ylabel('number of scattered photons')
xlabel(scan_var_name)
%  legendCell = cellstr(num2str(aaa',[regexprep(scan_var_name,'\',' ') '=%-f']))
%  legend(legendCell,'FontSize',16)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);




%% plot spectrum in the theta
nbin=50;
min_val_linspase=50
figure(ifig)
ifig=ifig+1;
hold on
cc=jet(ww);
for ni=1:1:ww
    plot(linspace(min_val_linspase,max(eval(['full_spectrum_' int2str(ni) ])),nbin),...
        histc(eval(['full_spectrum_' int2str(ni) '(aa_' int2str(ni) ')']),...
        linspace(min_val_linspase,max(eval(['full_spectrum_' int2str(ni) ])),nbin))*eval(['weigth_' int2str(ni) ]),'color',cc(ni,:),'LineWidth',2)
    
    bwt(ni)=eval(['std(full_spectrum_' int2str(ni) '(aa_' int2str(ni) '))/mean(full_spectrum_' int2str(ni) '(aa_' int2str(ni) '))'])
end
grid on
set(gca,'FontSize',16)
title({['Spectrum of scattered photons in Theta=' num2str(theta) ' [rad]'],...
    ['bandwith =' num2str(std(full_spectrum_1(aa_1))/mean(full_spectrum_1(aa_1))) ...
    '; ' num2str(std(full_spectrum_2(aa_2))/mean(full_spectrum_2(aa_2))) ';'...
    num2str(std(full_spectrum_3(aa_3))/mean(full_spectrum_3(aa_3))) '']})
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')
legendCell = cellstr([num2str(aaa', 'def=%-10.1f; ') num2str(bwt', ' BW=%-f')])
legend(legendCell,'Location','northwest','FontSize',16)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);



figure(ifig)
ifig=ifig+1;
subplot 211
plot(scan_var,bwt,'LineWidth',2)
set(gca,'FontSize',16)
ylabel('BW')
grid on
subplot 212
plot(scan_var,num_in_b,'LineWidth',2)
grid on
set(gca,'FontSize',16)
ylabel('Num in BW')
xlabel(scan_var_name)
suptitle(['For \theta=' num2str(theta) ])
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);
% title({['Spectrum of scattered photons in Theta=' num2str(theta) ' [rad]'],...
%      ['bandwith =' num2str(std(full_spectrum_1(aa_1))/mean(full_spectrum_1(aa_1))) ...
%      '; ' num2str(std(full_spectrum_2(aa_2))/mean(full_spectrum_2(aa_2))) ';'...
%       num2str(std(full_spectrum_3(aa_3))/mean(full_spectrum_3(aa_3))) '']})
% ylabel('number of scattered photons')
% xlabel('photons energy (KeV)')
%  legendCell = cellstr([num2str(aaa', 'def=%-10.1f; ') num2str(bwt', ' BW=%-f')])
%  legend(legendCell,'Location','northwest','FontSize',16)
% filename = ['plot_' num2str(ifig) ];
% fname = [ filename '.png'];
% print('-dpng', fname);



%% find theta for bandwidth 0.5% for each event

ww=0;
for ni=aaa
    ww=ww+1
    bandwith_cm=[];num_phot_th=[];aam=[];
    theta_angle=5e-6;
    diapason=(10:1:1001);
    el_angel=0;
    qq=0; num_phot_th=[]; bandwith_cm=[]; bandwith_non_norm=[]; aa_1=[];el_angel_1=[];
    for ni=diapason;
        qq=qq+1;
        el_angel_1=ni*theta_angle;
        eval(['aa_bw=find(abs(phot_angle_' int2str(ww) ')<el_angel_1);'])
        
        num_phot_th(qq)=length(full_spectrum(aa_1))*weigth;
        bandwith_cm(qq)=eval(['std(full_spectrum_' int2str(ww) '(aa_bw))/mean(full_spectrum_' int2str(ww) '(aa_bw))']);
        if bandwith_cm(qq)>=0.00495 && bandwith_cm(qq)<=0.101 % find angle for given bandwith
            el_angel=el_angel_1;
            % stop;
        end
    end
    
    %%%%
    thetam(ww)=el_angel;%19.25e-5;
    eval(['aam=find(abs(phot_angle_' int2str(ww) ')<thetam(ww));']);
    eval(['bandwith_in_thetam(ww)=std(full_spectrum_' int2str(ww) '(aam))/mean(full_spectrum_' int2str(ww) '(aam));'])
    eval(['num_in_thetam(ww)=length(find(full_spectrum_' int2str(ww) '(aam)))*weigth;'])
    eval(['BW_EV(ww)=std(1e3*full_spectrum_' int2str(ww) '(aam));'])
    
    eval(['brilliance_peak(ww)= num_in_thetam(ww)/((std(z_coor' int2str(ww) '(aam))/3e8)*((2*pi)^(5/2))*std(x_coor' int2str(ww) '(aam))*1e3*std(xp_' int2str(ww) '(aam)./zp_' int2str(ww) '(aam))*1e3*std(y_coor' int2str(ww) '(aam))*1e3*std(yp_' int2str(ww) '(aam)./zp_' int2str(ww) '(aam))*1e3*5);'])
    
%     eval(['colim_spectr(ww)=full_' int2str(ww) 'spectrum(aam);'])
%     eval(['in_the_sigma(ww)=find(colim_spectr(ww)>(pick_of_spectr-(1/2)*std(full_' int2str(ww) 'spectrum(aam))) & colim_spectr<(pick_of_spectr+(1/2)*std(full_' int2str(ww) 'spectrum(aam))));'])
    
      
end






figure(ifig)
ifig=ifig+1;
subplot 211
plot(scan_var,brilliance_peak,'-*b','LineWidth',2)
grid on
set(gca,'FontSize',16)
ylabel('brilliance peak')
subplot 212
plot(scan_var,3200*num_in_thetam./(sqrt(2*pi).*BW_EV),'-*g','LineWidth',2)
grid on
set(gca,'FontSize',16)
ylabel('SPD')
xlabel(scan_var_name)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-r300','-dpng', fname);


figure(ifig)
ifig=ifig+1;
subplot 411
plot(scan_var,bandwith_in_thetam,'-*b','LineWidth',2)
grid on
set(gca,'FontSize',16)
ylabel('BW')
% xlabel('defocus param')
subplot 412
plot(scan_var,num_in_thetam,'-*r','LineWidth',2)
grid on
set(gca,'FontSize',16)
ylabel('num in bw')
subplot 413
plot(scan_var,thetam,'-*g','LineWidth',2)
grid on
set(gca,'FontSize',16)
ylabel('\theta_{max}')
subplot 414
plot(scan_var,num_in_thetam./(sqrt(2*pi).*BW_EV),'-*g','LineWidth',2)
grid on
set(gca,'FontSize',16)
title('num in BW/(sqrt(2*pi).*BW(eV))')
ylabel('SPD')
xlabel(scan_var_name)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-r300','-dpng', fname);



figure(ifig)
ifig=ifig+1;
hold on
plot(scan_var,NUM_ph_x,'--xr','LineWidth',2)
plot(scan_var,total_num_in_b,'--ob','LineWidth',2)
hold off
 set(gca,'FontSize',22)
legend(' Analytical','Simulations cain',0)
xlabel([scan_var_name ' [deg]'])%,'FontSize',20)
ylabel('Number of the scattering photons')%,'FontSize',20)
grid on
% title(['\sigma_x=' num2str(std_x_size*1e6) ' \mu m; \sigma_y=' num2str(std_y_size*1e6) ' \mu m; \sigma_z=' num2str(std_s*1e3) ' mm;'])
filename = [ 'L_Flux_vs_alpha_' num2str(500)];
fname = [ filename '.png'];
print('-dpng', fname);
