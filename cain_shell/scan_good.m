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
for var_for_scan=[18 20 26 28 30 32]%0:10:100
    qq=qq+1;
    rflags.sigLr=var_for_scan/2; % given in [mu m] micro meter like 2 weist w0=28;
    
    defoc_param=1.6;%defocusing parameter by defoult 1 is no defocusing
    BASE_DIRECTORY = [pwd '/hi_400_defocus_' num2str(defoc_param) '_waist_' num2str(var_for_scan) '/'];
    
    [full_spectrum,phot_angle,weigth]=scan_photons_plots(BASE_DIRECTORY);
        
      eval(['full_spectrum_' int2str(qq) '=full_spectrum;'])
      eval(['phot_angle_' int2str(qq) '=phot_angle;'])
    eval(['weigth_' int2str(qq) '=weigth;'])
end


% color_line=['-r','-b','-g','-m','-y','-k'];
cc=jet(qq);

ifig=1;
figure(ifig)
ifig=ifig+1;
for ni=1:1:qq
hold on
plot(linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))),...
    histc(eval(['full_spectrum_' int2str(ni) ]),linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))))*eval(['weigth_' int2str(qq) ]),...
    'color',cc(ni,:),'LineWidth',2)
grid on
hold off
end
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')
legend('waist=26', 'waist=28', 'waist=30', 'waist=32',0)% give by hand sorry
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);



el_angel=0.0002;
aa_1=find(abs(phot_angle_1)<el_angel);
aa_2=find(abs(phot_angle_2)<el_angel);
aa_3=find(abs(phot_angle_3)<el_angel);
aa_4=find(abs(phot_angle_4)<el_angel);

nbin=50;
min_val_linspase=2600
figure(ifig)
ifig=ifig+1;
hold on
plot(linspace(min_val_linspase,max(full_spectrum_1),nbin),histc(full_spectrum_1(aa_1),linspace(min_val_linspase,max(full_spectrum_1),nbin))*weigth,'-r','LineWidth',2)
plot(linspace(min_val_linspase,max(full_spectrum_2),nbin),histc(full_spectrum_2(aa_2),linspace(min_val_linspase,max(full_spectrum_2),nbin))*weigth,'-b','LineWidth',2)
plot(linspace(min_val_linspase,max(full_spectrum_3),nbin),histc(full_spectrum_3(aa_3),linspace(min_val_linspase,max(full_spectrum_3),nbin))*weigth,'-g','LineWidth',2)
plot(linspace(min_val_linspase,max(full_spectrum_4),nbin),histc(full_spectrum_4(aa_4),linspace(min_val_linspase,max(full_spectrum_4),nbin))*weigth,'-m','LineWidth',2)
grid on
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')
legend('waist=26', 'waist=28', 'waist=30', 'waist=32',0)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);
