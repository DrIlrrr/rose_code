% one beam several events
clear all; close all; clc;
make_path
start_date=datestr(now);
%for new global
global home_dir DIRECTORY_FOR_CAIN BASE_DIRECTORY;
global rflags
[rflags] = flags_for_run;
home_dir=[pwd '/CAIN/'];

rflags.PLOTS =1;
just_plots=1;
qq=0;
scan_spectrum=[];
phot_angle=[];



aaa=[0.0001 0.2 0.4 1 4];


for var_for_scan=aaa
    qq=qq+1;
    
   defoc_param=1.5;%defocusing parameter by defoult 1 is no defocusing
   BASE_DIRECTORY = [pwd '/dif_polar_pulseE_' num2str(var_for_scan) '_hi_600_defocus_' num2str(defoc_param) '_PulseE_400mJ/'];


    [full_spectrum,phot_angle,weigth]=scan_photons_plots(BASE_DIRECTORY);
       
      eval(['full_spectrum_' int2str(qq) '=full_spectrum;'])
      eval(['phot_angle_' int2str(qq) '=phot_angle;'])
    
    
    
end


color_line={'-r','-b','-g','-m','-y'};



ifig=1;
figure(ifig)
ifig=ifig+1;
for ni=1:1:4
hold on
plot(linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))),...
    histc(eval(['full_spectrum_' int2str(ni) ]),linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))))*weigth,...
    color_line{ni},'LineWidth',2)
% plot(linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))),histc(eval(['full_spectrum_' int2str(ni) ]),linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))))*weigth,color_line(ni),'LineWidth',2)
% plot(linspace(0,max(full_spectrum_1)),histc(full_spectrum_1,linspace(0,max(full_spectrum_1)))*weigth,'-r','LineWidth',2)
grid on
hold off
end
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')
legend('\delta t=0 [ps]','\delta t=1 [ps]','\delta t=2 [ps]',0)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);
% 
% ifig=1;
% figure(ifig)
% ifig=ifig+1;
% hold on
% plot(linspace(0,max(full_spectrum_1)),histc(full_spectrum_1,linspace(0,max(full_spectrum_1)))*weigth,'-r','LineWidth',2)
% plot(linspace(0,max(full_spectrum_2)),histc(full_spectrum_2,linspace(0,max(full_spectrum_2)))*weigth,'-b','LineWidth',2)
% plot(linspace(0,max(full_spectrum_3)),histc(full_spectrum_3,linspace(0,max(full_spectrum_3)))*weigth,'-g','LineWidth',2)
% plot(linspace(0,max(full_spectrum_4)),histc(full_spectrum_4,linspace(0,max(full_spectrum_4)))*weigth,'-m','LineWidth',2)
% grid on
% set(gca,'FontSize',16)
% ylabel('number of scattered photons')
% xlabel('photons energy (KeV)')
% legend('/del=26', 'waist=28', 'waist=30', 'waist=32',0)
% filename = ['plot_' num2str(ifig) ];
% fname = [ filename '.png'];
% print('-dpng', fname);


el_angel=8.8e-5;
% aa_1=find(abs(phot_angle_1)<el_angel);
% aa_2=find(abs(phot_angle_2)<el_angel);
% aa_3=find(abs(phot_angle_3)<el_angel);
% aa_4=find(abs(phot_angle_4)<el_angel);
% aa_5=find(abs(phot_angle_5)<el_angel);
% aa_6=find(abs(phot_angle_6)<el_angel);
% aa_7=find(abs(phot_angle_7)<el_angel);

aa_1=find(abs(phot_angle_1));
aa_2=find(abs(phot_angle_2));
aa_3=find(abs(phot_angle_3));
aa_4=find(abs(phot_angle_4));
aa_5=find(abs(phot_angle_5));
aa_6=find(abs(phot_angle_6));
aa_7=find(abs(phot_angle_7));




ww=0;
for ni=aaa
    ww=ww+1
    num_in_b(ww)=eval(['length(full_spectrum_' int2str(ww) '(aa_' int2str(ww) '))']);
    
end

figure(ifig)
ifig=ifig+1;
hold on
plot(aaa,num_in_b*weigth,'-.xb','LineWidth',0.5)
plot(aaa,num_in_b*weigth,'xb','LineWidth',3)
hold off
grid on
ylim([0 max(num_in_b*weigth)+max(num_in_b*weigth)*1e-1])
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('time delay [ps]')
% legend('\delta t=0 [ps]','\delta t=1 [ps]','\delta t=2 [ps]',0)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);



nbin=50;
min_val_linspase=2600
figure(ifig)
ifig=ifig+1;
hold on
plot(linspace(min_val_linspase,max(full_spectrum_1),nbin),histc(full_spectrum_1(aa_1),linspace(min_val_linspase,max(full_spectrum_1),nbin))*weigth,'-k','LineWidth',2)
plot(linspace(min_val_linspase,max(full_spectrum_2),nbin),histc(full_spectrum_2(aa_2),linspace(min_val_linspase,max(full_spectrum_2),nbin))*weigth,'--b','LineWidth',2)
plot(linspace(min_val_linspase,max(full_spectrum_3),nbin),histc(full_spectrum_3(aa_3),linspace(min_val_linspase,max(full_spectrum_3),nbin))*weigth,'-g','LineWidth',2)
plot(linspace(min_val_linspase,max(full_spectrum_4),nbin),histc(full_spectrum_4(aa_4),linspace(min_val_linspase,max(full_spectrum_4),nbin))*weigth,'-m','LineWidth',2)
plot(linspace(min_val_linspase,max(full_spectrum_5),nbin),histc(full_spectrum_5(aa_5),linspace(min_val_linspase,max(full_spectrum_5),nbin))*weigth,'-c','LineWidth',2)
plot(linspace(min_val_linspase,max(full_spectrum_6),nbin),histc(full_spectrum_6(aa_6),linspace(min_val_linspase,max(full_spectrum_6),nbin))*weigth,'-y','LineWidth',2)
plot(linspace(min_val_linspase,max(full_spectrum_7),nbin),histc(full_spectrum_7(aa_7),linspace(min_val_linspase,max(full_spectrum_7),nbin))*weigth,'-r','LineWidth',2)
grid on
set(gca,'FontSize',16)
title({['Spectrum of scattered photons in Theta=' num2str(el_angel) ' [rad]'],...
    ['bandwith =' num2str(std(full_spectrum_1(aa_1))/mean(full_spectrum_1(aa_1))) ' for \delta t=0']})
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')
legend(['\delta t=0 [ps]; photons = ' num2str(num_in_b(1)*weigth,'%2.1e') ],...
    ['\delta t=4 [ps]; photons = ' num2str(num_in_b(2)*weigth,'%2.1e') ],...
    ['\delta t=5 [ps]; photons = ' num2str(num_in_b(3)*weigth,'%2.1e') ],...
   ['\delta t=6 [ps]; photons = ' num2str(num_in_b(4)*weigth,'%2.1e') ],...
   ['\delta t=7 [ps]; photons = ' num2str(num_in_b(5)*weigth,'%2.1e') ],...
   ['\delta t=8 [ps]; photons = ' num2str(num_in_b(6)*weigth,'%2.1e') ],...
   ['\delta t=9 [ps]; photons = ' num2str(num_in_b(7)*weigth,'%2.1e') ],0)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);





% figure(ifig)
% ifig=ifig+1;
% for ni=1:1:3
% hold on
% plot((eval(['full_spectrum_' int2str(1:1:3) ])(aa_1))*weigth,...
%     color_line{ni},'LineWidth',2)
% % plot(linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))),histc(eval(['full_spectrum_' int2str(ni) ]),linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))))*weigth,color_line(ni),'LineWidth',2)
% % plot(linspace(0,max(full_spectrum_1)),histc(full_spectrum_1,linspace(0,max(full_spectrum_1)))*weigth,'-r','LineWidth',2)
% grid on
% hold off
% end
% set(gca,'FontSize',16)
% ylabel('number of scattered photons')
% xlabel('photons energy (KeV)')
% legend('\delta t=0 [ps]','\delta t=1 [ps]','\delta t=2 [ps]',0)
% filename = ['plot_' num2str(ifig) ];
% fname = [ filename '.png'];
% print('-dpng', fname);


