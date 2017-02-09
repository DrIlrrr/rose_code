clear all; close all; clc;
make_path
start_date=datestr(now);
%for new global
% global home_dir DIRECTORY_FOR_CAIN BASE_DIRECTORY;
% global rflags
[rflags] = flags_for_run;
home_dir=[pwd '/CAIN/'];
mkdir(home_dir)
diapason=[0];
index_num=0;


for aa_sf=diapason
    index_num=index_num+1;
rflags.shifting_laser_t=1e-12*3e8*aa_sf;
close all;
% out_folder = [pwd '/try_long_laser_shifting_' num2str(rflags.shifting_laser_t) '/'];
out_folder = [pwd '/new_beam/'];

electron_data=load([out_folder 'cain_tmp/exp.dat']);
std_x_size=std(electron_data(:,5));
std_y_size=std(electron_data(:,6));
std_s=std(electron_data(:,4));
number_electrons=electron_data(1,3)*length(electron_data(:,3));

SPEED_OF_LIGHT=3e8;
Sigma_th=0.665e-28;%[m^2] Thomsom cross section
h=2*pi*1.054e-34;%Planc const [J*s]
lambda_l=515e-9;% [m] laser wave length
pulseE=0.2; %[J]
photons_number=pulseE/((h*SPEED_OF_LIGHT)/lambda_l);

angle=7.5;
sigLr=28e-6/2.0;
%         laser_length=var_val;
sigt=1.5e-12;
laser_length=sigt*SPEED_OF_LIGHT;
delta_s=1e-12*3e8*aa_sf;
delta_x=0;

As=exp((delta_s^2.*tan((1/2)*angle.*(pi/180)))/(-2*((std_x_size.^2+sigLr^2)+((std_s).^2+(sigt.*SPEED_OF_LIGHT).^2).*tan((1/2)*angle.*(pi/180)))));
Ax=exp((delta_x^2)/(-2*((std_x_size.^2+sigLr^2)+((std_s).^2+(sigt.*SPEED_OF_LIGHT).^2).*tan((1/2)*angle.*(pi/180)))));




NUM_ph(index_num)=As.*Ax.*((Sigma_th*photons_number*number_electrons)/(2*pi*sqrt(std_y_size.^2+sigLr^2)))...
    /sqrt((std_x_size.^2+sigLr^2)+((std_s).^2+(sigt.*SPEED_OF_LIGHT).^2).*tan((1/2)*angle.*(pi/180))^2);

for turn_number=[3]
mkdir([ out_folder 'photon_plots_' num2str(turn_number)])

%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss

all_data_phot=[];
full_spectrum=[];
phot_angle=[];
weigth=0;
histmat=[];
new_m=[];
x_phot=[];
y_phot=[];
xp_phot=[];
yp_phot=[];

for ni=1:1:turn_number
   
    
    
    
check_photons_output=dir([out_folder 'cain_tmp/cain_output_photons_' num2str(ni) '.dat']);% checking cain output for photons

check_photons_output_size=check_photons_output(:).bytes;

if (check_photons_output_size==0)%if cain non produce photons write photons_property=[0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    photons_data=zeros(1,14);
else
    
 photons_data=dlmread([out_folder 'cain_tmp/cain_output_photons_' num2str(ni) '.dat'],'',1,0);%read electrons from cain
end
    
   
    
    x_phot=[x_phot;photons_data(:,5)];
    y_phot=[y_phot;photons_data(:,6)];
    xp_phot=[xp_phot;photons_data(:,9)];
    yp_phot=[yp_phot;photons_data(:,10)];
    
    weigth=photons_data(1,3);
%     phot_number=size(photons_data,1)*weigth;
%     zx=size(photons_data(:,8));
%     number_of_photons=zx(1)*weigth;
    
    
    full_spectrum=[full_spectrum;photons_data(:,8)./1e3];
    % phot_angle=[phot_angle;sign(photons_data(:,9)).*atan(sqrt(photons_data(:,9).^2+photons_data(:,10).^2)./photons_data(:,11))];
    phot_angle=[phot_angle;atan(sqrt(photons_data(:,9).^2+photons_data(:,10).^2)./photons_data(:,11))];
    
end

number_of_photons=length(full_spectrum)*weigth;
%SPEED_OF_LIGHT=3e8;
el_angel=4e-5;
aa=find(abs(phot_angle)<el_angel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
number_of_photons_shifted(index_num)=length(full_spectrum)/turn_number*weigth;
%SPEED_OF_LIGHT=3e8;
el_angel_m=4e-5;
aa=find(abs(phot_angle)<el_angel_m);
number_of_photons_shifted_in_angle(index_num)=length(full_spectrum(aa))/turn_number*weigth;
bandwith_cm_shifting(index_num)=std(full_spectrum(aa))/mean(full_spectrum(aa));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ifig=1;

bandwith_cm=[];
num_phot_th=[];
for ni=1:1:100
el_angel=ni*0.5e-6;
aa_1=find(abs(phot_angle)<el_angel);

num_phot_th(ni)=length(full_spectrum(aa_1))/turn_number*weigth;
bandwith_cm(ni)=std(full_spectrum(aa_1))/mean(full_spectrum(aa_1));

end
no_plots=1;
if no_plots==0
figure(ifig)
ifig=ifig+1;
subplot(2,1,1)
plot((1:1:100).*0.5e-6,num_phot_th,'-xb')
grid on
ylabel('number photons ')
subplot(2,1,2)
plot((1:1:100).*0.5e-6,bandwith_cm,'--sb')
ylabel('bandwith')
grid on
filename = [ out_folder 'photon_plots_' num2str(turn_number) '/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);


figure(ifig)
ifig=ifig+1;
subplot(2,1,1)
plot((1:1:100).*0.5e-6,num_phot_th,'-xb')
grid on
ylabel('number photons ')
subplot(2,1,2)
plot((1:1:100).*0.5e-6,bandwith_cm,'--sb')
ylabel('bandwith')
grid on
suptitle([ 'Number events ' num2str(turn_number) ])
filename = [ out_folder 'fig_4_num_ev_' num2str(turn_number) ];
fname = [ filename '.png'];
print('-dpng', fname);




figure(ifig)
ifig=ifig+1;
subplot(1,2,1)
hold on
hist(phot_angle,50)
hold off
grid on
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('theta [rad]')
% xlabel('photons energy (KeV)')
subplot(1,2,2)
hold on
hist(full_spectrum,50)
hold off
grid on
set(gca,'FontSize',16)
ylabel('number of scattered photons')
%xlabel('theta')
 xlabel('photons energy (KeV)')
suptitle({['For all theta'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ]})
filename = [ out_folder 'photon_plots_' num2str(turn_number) '/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);



figure(ifig)
ifig=ifig+1;
subplot(1,2,1)
hold on
hist(phot_angle(aa),50)
hold off
grid on
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('theta')
% xlabel('photons energy (KeV)')
subplot(1,2,2)
hold on
hist(full_spectrum(aa),50)
hold off
grid on
set(gca,'FontSize',16)
ylabel('number of scattered photons')
%xlabel('theta')
 xlabel('photons energy (KeV)')
 suptitle({['For theta<' num2str(el_angel) ' [rad]'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ];['number scattered photons in theta<' num2str(length(full_spectrum(aa))*weigth/turn_number,'%10.2e') ]})

filename = [ out_folder 'photon_plots_' num2str(turn_number) '/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);



nbin_plot=30;
figure(ifig)
ifig=ifig+1;
hold on
plot(linspace(0,max(full_spectrum),nbin_plot),smooth(hist(full_spectrum,nbin_plot)*weigth),'-r','LineWidth',2)
hold off
grid on
ylim([0 max(smooth(hist(full_spectrum,nbin_plot)*weigth))])
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')

title({['For all theta'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ]})
filename = [ out_folder 'photon_plots_' num2str(turn_number) '/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);




figure(ifig)
ifig=ifig+1;
hold on
plot(phot_angle,full_spectrum,'.b')
hold off
grid on
set(gca,'FontSize',16)
ylabel('energy of scattered photons')
xlabel('scattered angle')
suptitle({['For all theta'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ]});
filename = [ out_folder 'photon_plots_' num2str(turn_number) '/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);






figure(ifig)
ifig=ifig+1;
hold on
plot(phot_angle(aa),full_spectrum(aa),'.b')
hold off
grid on
set(gca,'FontSize',16)
ylabel('energy of scattered photons')
xlabel('scattered angle')
suptitle({['For theta<' num2str(el_angel) ' [rad]'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ];['number scattered photons in theta<' num2str(length(full_spectrum(aa))*weigth/turn_number,'%10.2e') ]})

filename = [ out_folder 'photon_plots_' num2str(turn_number) '/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);

figure(ifig)
ifig=ifig+1;
hold on
plot(linspace(0,max(full_spectrum)),histc(full_spectrum(aa),linspace(0,max(full_spectrum)))*weigth,'-r','LineWidth',2)
% plot([30:0.5:47],histc(full_spectrum(aa),[30:0.5:47])*weigth,'-r','LineWidth',2)
hold off
grid on
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')
% suptitle(['For theta<' num2str(el_angel) ' [rad]']);
suptitle({['For theta<' num2str(el_angel) ' [rad]'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ];['number scattered photons in theta<' num2str(length(full_spectrum(aa))*weigth/turn_number,'%10.2e') ]})
text(93.1658986175114, 110886076.949367,{['bandwidth=' num2str(std(full_spectrum(aa))/mean(full_spectrum(aa)))];...
    ['FWHM=' num2str(0)]})
filename = [ out_folder 'photon_plots_' num2str(turn_number) '/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);









ifig=ifig+1;
figure(ifig)
subplot 221
set(gca,'FontSize',16)
plot(x_phot,y_phot,'.b')
xlabel('x [m]')
ylabel('y [m]')
subplot 222
set(gca,'FontSize',16)
plot(x_phot,xp_phot,'.b')
xlabel('x [m]')
ylabel('Px [eV/c]')
subplot 223
set(gca,'FontSize',16)
plot(y_phot,yp_phot,'.b')
xlabel('y [m]')
ylabel('Px [eV/c]')
subplot 224
set(gca,'FontSize',16)
plot(xp_phot,yp_phot,'.b')
xlabel('Px [eV/c]')
ylabel('Py [eV/c]')
suptitle({['For all theta'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ]});
filename = [ out_folder 'photon_plots_' num2str(turn_number) '/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);



ifig=ifig+1;
figure(ifig)

subplot 221
set(gca,'FontSize',16)
plot(x_phot(aa),y_phot(aa),'.b')
xlabel('x [m]')
ylabel('y [m]')
subplot 222
set(gca,'FontSize',16)
plot(x_phot(aa),xp_phot(aa),'.b')
xlabel('x [m]')
ylabel('Px [eV/c]')
subplot 223
set(gca,'FontSize',16)
plot(y_phot(aa),yp_phot(aa),'.b')
xlabel('y [m]')
ylabel('Px [eV/c]')
subplot 224
set(gca,'FontSize',16)
plot(xp_phot(aa),yp_phot(aa),'.b')
xlabel('Px [eV/c]')
ylabel('Py [eV/c]')

suptitle({['For theta<' num2str(el_angel) ' [rad]'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ];['number scattered photons in theta<' num2str(length(full_spectrum(aa))*weigth/turn_number,'%10.2e') ]})


filename = [ out_folder 'photon_plots_' num2str(turn_number) '/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);








ifig=ifig+1;
figure(ifig)
subplot 221

set(gca,'FontSize',16)
hold on
plot(x_phot,y_phot,'.b')
plot(x_phot(aa),y_phot(aa),'.r')
hold off
xlabel('x [m]')
ylabel('y [m]')
subplot 222
hold on
set(gca,'FontSize',16)
plot(x_phot,xp_phot,'.b')
plot(x_phot(aa),xp_phot(aa),'.r')
hold off
xlabel('x [m]')
ylabel('Px [eV/c]')
subplot 223
hold on
set(gca,'FontSize',16)
plot(y_phot,yp_phot,'.b')
plot(y_phot(aa),yp_phot(aa),'.r')
hold off
xlabel('y [m]')
ylabel('Px [eV/c]')
subplot 224
set(gca,'FontSize',16)
hold on
plot(xp_phot,yp_phot,'.b')
plot(xp_phot(aa),yp_phot(aa),'.r')

hold off
xlabel('Px [eV/c]')
ylabel('Py [eV/c]')

suptitle({['For theta<' num2str(el_angel) ' [rad]'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ];['number scattered photons in theta<' num2str(length(full_spectrum(aa))*weigth/turn_number,'%10.2e') ]})

filename = [ out_folder 'photon_plots_' num2str(turn_number) '/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);

end

end

end
% 
% number_of_photons_shifted(aa_sf)=length(full_spectrum)/turn_number*weigth;
% %SPEED_OF_LIGHT=3e8;
% el_angel=4e-5
% aa=find(abs(phot_angle)<el_angel);
% number_of_photons_shifted_in_angle(aa_sf)=length(full_spectrum(aa))/turn_number*weigth;
% bandwith_cm(aa_sf)=std(full_spectrum(aa))/mean(full_spectrum(aa));



ifig=ifig+1;
figure(ifig)
subplot 311
set(gca,'FontSize',14)
hold on
 plot(1e-12*3e8*diapason,number_of_photons_shifted,'-ob')
  plot(1e-12*3e8*diapason,NUM_ph,'-or')
 hold off
 for nii=1:1:length(diapason)
 text(1e-12*3e8*diapason(nii),number_of_photons_shifted(nii),['' num2str(number_of_photons_shifted(nii),'%10.0e\')],'FontSize',14) 
 text(1e-12*3e8*diapason(nii),NUM_ph(nii),['' num2str(NUM_ph(nii),'%10.0e\')],'FontSize',14) 
 end
 grid on
%plot(number_of_photons_shifted,'.b')
xlabel('s [m]')
ylabel('Full number of photons')
subplot 312
set(gca,'FontSize',14)
 plot(1e-12*3e8*diapason,number_of_photons_shifted_in_angle,'-ob')
% plot(number_of_photons_shifted_in_angle,'.b')
 grid on
 for nii=1:1:length(diapason)
 text(1e-12*3e8*diapason(nii),number_of_photons_shifted_in_angle(nii),['' num2str(number_of_photons_shifted_in_angle(nii),'%10.0e\')],'FontSize',14) 
 end
xlabel('s [m]')
ylabel(['Number of photons in theta=' num2str(el_angel_m)])
subplot 313
set(gca,'FontSize',14)
plot(1e-12*3e8*diapason,bandwith_cm_shifting,'-ob')
% plot(bandwith_cm_shifting,'.b')
 grid on
xlabel('s [m]')
ylabel('bandwith')
% suptitle({['For all theta'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ]});
 filename = ['photons_plot_' num2str(ifig) ];
 fname = [ filename '.png'];
 print('-dpng', fname);

