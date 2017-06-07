% one beam several events
clear all; close all; clc;
make_path
start_date=datestr(now);
%for new global
global home_dir DIRECTORY_FOR_CAIN BASE_DIRECTORY;
global rflags
[rflags] = flags_for_run;
home_dir=[pwd '/CAIN/'];
ifig=1;
rflags.PLOTS =1;
just_plots=1;
qq=0;
scan_spectrum=[];
phot_angle=[];


aaa=[0:2.5:7.5]


for var_for_scan=aaa
    qq=qq+1;
    
   defoc_param=1.5;%defocusing parameter by defoult 1 is no defocusing
   BASE_DIRECTORY = [pwd '/dif_angle_' num2str(var_for_scan) '_hi_600_defocus_' num2str(defoc_param) '_PulseE_400mJ/'];
    [full_spectrum,phot_angle,weigth]=scan_photons_plots(BASE_DIRECTORY);
       
      eval(['full_spectrum_' int2str(qq) '=full_spectrum;'])
      eval(['phot_angle_' int2str(qq) '=phot_angle;'])
    
   electron_data=load([BASE_DIRECTORY 'cain_tmp/exp.dat']);
    
   std_x_size=std(electron_data(:,5));
std_y_size=std(electron_data(:,6));
std_s=std(electron_data(:,4));
number_electrons=electron_data(1,3)*length(electron_data(:,3));

SPEED_OF_LIGHT=3e8;
Sigma_th=0.665e-28;%[m^2] Thomsom cross section
h=2*pi*1.054e-34;%Planc const [J*s]
lambda_l=515e-9;% [m] laser wave length
pulseE=0.4; %[J]
photons_number=pulseE/((h*SPEED_OF_LIGHT)/lambda_l);

angle=var_for_scan;
sigLr=28e-6/2.0;
%         laser_length=var_val;
sigt=1.5e-12;
laser_length=sigt*SPEED_OF_LIGHT;


% As=exp((delta_s^2.*tan((1/2)*angle.*(pi/180)))/(-2*((std_x_size.^2+sigLr^2)+((std_s).^2+(sigt.*SPEED_OF_LIGHT).^2).*tan((1/2)*angle.*(pi/180)))));
% Ax=exp((delta_x^2)/(-2*((std_x_size.^2+sigLr^2)+((std_s).^2+(sigt.*SPEED_OF_LIGHT).^2).*tan((1/2)*angle.*(pi/180)))));
%
%

% 
NUM_ph(qq)=((Sigma_th*photons_number*number_electrons)/(2*pi*sqrt(std_y_size^2+sigLr^2)))...
    /sqrt((std_x_size^2+sigLr^2)+((std_s)^2+(sigt*SPEED_OF_LIGHT)^2)*tan((1/2)*angle*(pi/180))^2);
   
   




bandwith_cm=[];
num_phot_th=[];
theta_angle=1e-6;
diapason=(10:1:101);
wqq=0; num_phot_th=[]; bandwith_cm=[]; bandwith_non_norm=[];
for ni=diapason;
    wqq=wqq+1;
    el_angel_1=ni*theta_angle;
    aa_1=find(abs(phot_angle)<el_angel_1);
    
    num_phot_th(wqq)=length(full_spectrum(aa_1))*weigth;
    bandwith_cm(wqq)=std(full_spectrum(aa_1))/mean(full_spectrum(aa_1));
    bandwith_non_norm(wqq)=std(full_spectrum(aa_1));
    
end

figure(ifig)
ifig=ifig+1;
set(gca,'FontSize',16)
subplot(2,1,1)
plot(diapason.*theta_angle,num_phot_th,'-xb')
grid on
ylabel('number photons ')
subplot(2,1,2)
plot(diapason.*theta_angle,bandwith_cm,'--sb')
ylabel('bandwith')
xlabel('Theta')
grid on
suptitle(['\alpha_0=' num2str(var_for_scan) ])
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);









   
end


color_line={'-r','-b','-g','-m','-y'};




figure(ifig)
ifig=ifig+1;
for ni=1:1:qq
hold on
plot(linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))),...
    histc(eval(['full_spectrum_' int2str(ni) ]),linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))))*weigth,...
    color_line{ni},'LineWidth',2)

% plot(linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))),...
%     histc(eval(['full_spectrum_' int2str(ni) ]),linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))))*weigth,...
%     'LineWidth',2);

grid on
hold off
end
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')
legend('\alpha_0=0','\alpha_0=2.5','\alpha_0=5','\alpha_0=7.5',0)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);


ww=0;
for ni=aaa
    ww=ww+1
    eval(['aa_' num2str(ww) '=find(abs(phot_angle_' num2str(ww) '));']);
    num_in_b(ww)=eval(['length(full_spectrum_' int2str(ww) '(aa_' int2str(ww) '));']);
    
end

figure(ifig)
ifig=ifig+1;
hold on
plot(aaa,num_in_b*weigth,'-.xb','LineWidth',0.5)
plot(aaa,num_in_b*weigth,'xb','LineWidth',3)
plot(aaa,NUM_ph,'or','LineWidth',3)
hold off
grid on
ylim([0 max(num_in_b*weigth)+max(num_in_b*weigth)*1e-1])
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('\alpha_0')
% legend('\delta t=0 [ps]','\delta t=1 [ps]','\delta t=2 [ps]',0)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);







ww=0;
for ni=aaa
    ww=ww+1
    eval(['aa_' num2str(ww) '=find(abs(phot_angle_' num2str(ww) ')<7.8e-5);']);
    num_in_b(ww)=eval(['length(full_spectrum_' int2str(ww) '(aa_' int2str(ww) '));']);
    
end

nbin=50;
min_val_linspase=12000
figure(ifig)
ifig=ifig+1;
hold on
cc=jet(ww)
for ni=1:1:ww
plot(linspace(min_val_linspase,max(eval(['full_spectrum_' int2str(ni) ])),nbin),...
    histc(eval(['full_spectrum_' int2str(ni) '(aa_' int2str(ni) ')']),...
    linspace(min_val_linspase,max(eval(['full_spectrum_' int2str(ni) ])),nbin))*weigth,'color',cc(ni,:),'LineWidth',2)
end
grid on
set(gca,'FontSize',16)
% title({['Spectrum of scattered photons in Theta=' num2str(el_angel) ' [rad]'],...
%     ['bandwith =' num2str(std(full_spectrum_1(aa_1))/mean(full_spectrum_1(aa_1))) ' for \delta t=0']})
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')
legend(['\alpha_0=0 ; photons = ' num2str(num_in_b(1)*weigth,'%2.1e') ],...
    ['\alpha_0=2.5; photons = ' num2str(num_in_b(2)*weigth,'%2.1e') ],...
    ['\alpha_0=5; photons = ' num2str(num_in_b(3)*weigth,'%2.1e') ],...
    ['\alpha_0=7.5; photons = ' num2str(num_in_b(4)*weigth,'%2.1e') ],0)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);














