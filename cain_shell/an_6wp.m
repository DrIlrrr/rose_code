%not a function
clear all; close all; clc;

for WP=4;


if WP==1
    wppath='2_7GeV';
elseif WP==2
    wppath='2_7GeV_no_energy_spread';
elseif WP==3
     wppath='2_7GeV_no_em';   
elseif WP==4
    wppath='2_7GeV_no_energy_spread_no_em';
elseif WP==5
    wppath='ideal_el_beam_3';
elseif WP==6
    wppath='ideal_el_beam_600';
end


dir0=pwd;
dir1=[dir0 '/' wppath];
addpath(dir1);%,dir3);


%  make_path

global rflags

[rflags] = flags_for_run;
ifig=1;
WEIGHT=3;
TIME_COORDINATE=4;%now we start use s(m) for caine
X_COORDINATE=5;
Y_COORDINATE=6;
ENERGY_OF_PARTICLE=8;
X_MOMENTUM=9;
Y_MOMENTUM=10;
S_MOMENTUM=11;
% POLARISATION: 12 13 14
N_COMPTON_HIT=15;
TURN_LAST_COMPTON_HIT=16;
%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss

photons_data=load([ wppath '/photon_data_plots/photons_data.dat']);%photons_data=load([out_folder 'photon_data_plots/photons_data.dat']);
beam_property=load([ wppath '/cain_tmp/exp.dat']);%beam_property=load([out_folder 'cain_tmp/exp.dat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weigth=photons_data(1,3);
full_spectrum=photons_data(:,8)./1e3;
phot_angle=sign(photons_data(:,9)).*atan(sqrt(photons_data(:,9).^2+photons_data(:,10).^2)./photons_data(:,11));

bandwith_cm=[];
num_phot_th=[];
theta_angle=10e-7;
diapason=(0:1:51);
el_angel=0;
qq=0; num_phot_th=[]; bandwith_cm=[]; bandwith_non_norm=[];
for ni=diapason;
    qq=qq+1;
    el_angel_1=ni*theta_angle;
    aa_1=find(abs(phot_angle)<el_angel_1);
    
    num_phot_th(qq)=length(full_spectrum(aa_1))*weigth;
    bandwith_cm(qq)=std(full_spectrum(aa_1))/mean(full_spectrum(aa_1));
    theta_max(qq)=ni*theta_angle;
    if bandwith_cm(qq)>=0.0048 && bandwith_cm(qq)<=0.0051 % find angle for given bandwith
        el_angel=el_angel_1;
    end
    
    
    
end

aa=find(abs(phot_angle)<el_angel);
bandwith=std(full_spectrum(aa))/mean(full_spectrum(aa));
el_angel;

figure(ifig)
ifig=ifig+1;
set(gca,'FontSize',16)
subplot(2,1,1)
plot(diapason.*theta_angle,num_phot_th,'-xb')
grid on
set(gca,'FontSize',16)
ylabel('number photons ')
subplot(2,1,2)
hold on
plot(diapason.*theta_angle,bandwith_cm,'--xb')
plot(el_angel,bandwith,'or','LineWidth',2)
hold off
set(gca,'FontSize',16)
ylabel('bandwidth')
xlabel('Theta')
grid on
set(gca,'FontSize',16)
suptitle({['bandwith=' num2str(bandwith) ' for \theta=' num2str(el_angel) ];...
    ['full flux=' num2str(length(find(full_spectrum))) ';'...
    ' Flux in \theta =' num2str(length(find(full_spectrum(aa))),'%10.2e') ]})

ph_x=(photons_data(:,X_COORDINATE));
ph_xp=(photons_data(:,X_MOMENTUM)./(photons_data(:,S_MOMENTUM)));
ph_y=(photons_data(:,Y_COORDINATE));
ph_yp=(photons_data(:,Y_MOMENTUM)./(photons_data(:,S_MOMENTUM)));
% phasespace(5,:)=photons_data(:,TIME_COORDINATE)';
% phasespace(6,:)=((photons_data(:,ENERGY_OF_PARTICLE)-mean(photons_data(:,ENERGY_OF_PARTICLE)))/mean(photons_data(:,ENERGY_OF_PARTICLE)))';


ph_emx_sqrt=sqrt(mean(ph_x.^2)*mean(ph_xp.^2)-mean(ph_x.*ph_xp)^2);
ph_emy_sqrt=sqrt(mean(ph_y.^2)*mean(ph_yp.^2)-mean(ph_y.*ph_yp)^2);

ph_emx_sqrt_in_bandw=sqrt(mean(ph_x(aa).^2)*mean(ph_xp(aa).^2)-mean(ph_x(aa).*ph_xp(aa))^2);
ph_emy_sqrt_in_bandw=sqrt(mean(ph_y(aa).^2)*mean(ph_yp(aa).^2)-mean(ph_y(aa).*ph_yp(aa))^2);

e_x=(beam_property(:,X_COORDINATE));
e_xp=(beam_property(:,X_MOMENTUM)./(beam_property(:,S_MOMENTUM)));
e_y=(beam_property(:,Y_COORDINATE));
e_yp=(beam_property(:,Y_MOMENTUM)./(beam_property(:,S_MOMENTUM)));
e_s=beam_property(:,TIME_COORDINATE)';
e_ens=((beam_property(:,ENERGY_OF_PARTICLE)-mean(beam_property(:,ENERGY_OF_PARTICLE)))/mean(beam_property(:,ENERGY_OF_PARTICLE)))';
e_E=beam_property(:,ENERGY_OF_PARTICLE);

e_emx_sqrt=sqrt(mean(e_x.^2)*mean(e_xp.^2)-mean(e_x.*e_xp)^2);
e_emy_sqrt=sqrt(mean(e_y.^2)*mean(e_yp.^2)-mean(e_y.*e_yp)^2);

sigma_0=std(e_x);
sigma_0y=std(e_y);
W0=rflags.sigLr*2*1e-6;
gamma=mean(e_E)/(0.511e6);
delta_gamma=std(e_E./(0.511e6));
norm_em_x=sqrt(gamma^2-1)*e_emx_sqrt;
norm_em_y=sqrt(gamma^2-1)*e_emy_sqrt;
laser_bandwidth=1e-3;
M2=1.2;
laserwl=rflags.laserwl*1e-9;
a_0=4.3*(laserwl/W0)*sqrt(rflags.pulseE/rflags.sigt);


 save(['data_WP_' num2str(WP) '.dat'],'theta_max','delta_gamma','gamma','norm_em_x','sigma_0',...
     'norm_em_y','sigma_0y','laser_bandwidth','M2','laserwl','W0','a_0',...
     'theta_max','bandwith_cm','e_emx_sqrt','e_emy_sqrt','e_E','full_spectrum','WP')

%just_formula



% %
% 
figure(ifig)
ifig=ifig+1;
% subplot(1,2,1)
set(gca,'FontSize',16)
hold on
plot(theta_max.^2.*gamma.^2,bandwith_cm,'--xb')
plot(theta_max.^2.*gamma.^2,theta_max.^2.*gamma.^2.*(1/sqrt(12))+bandwith_cm(1)+0.5e-3,'-r')
hold off
grid on
xlabel('\theta_{max}^2*\gamma^2')
ylabel('bandwith')
legend('Simulation','(\theta_{max}^2*\gamma^2)/sqrt(12)','Location','best')
% subplot(1,2,2)
% set(gca,'FontSize',16)
% hold on
 plot(theta_max,bandwith_cm./(theta_max.^2.*gamma.^2),'--xb')
 plot(theta_max,1/sqrt(12),'.r')
% hold off
% grid on
% xlabel('Theta max')
% ylabel('bandwith/(\theta_{max}^2*\gamma^2)')
 legend('Simulation','1/sqrt(12)','Location','best')
 filename = ['comp_' num2str(ifig) ];
 fname = [ filename '.png'];
 print('-dpng', fname);
 
 %stop

% just_formula

f_name=['WP_' num2str(WP) '.txt']
%  filename = sprintf('NEW_plot.txt',tilted_name);
fileID = fopen(f_name,'w');
fprintf(fileID,'Photon Emittances full beam:\n');
fprintf(fileID,'PH Emit X = %10.5e\n',ph_emx_sqrt);
fprintf(fileID,'PH Emit Y = %10.5e\n',ph_emy_sqrt);
fprintf(fileID,'Photon Emittances collimated in theta = %10.5e with bandwidth = %10.5e:\n',el_angel,bandwith);
fprintf(fileID,'PH Emit X = %10.5e\n',ph_emx_sqrt_in_bandw);
fprintf(fileID,'PH Emit Y = %10.5e\n',ph_emy_sqrt_in_bandw);
fprintf(fileID,'----------------------------------------------------\n');
fprintf(fileID,'Electron beam non norm Emittances:\n');
fprintf(fileID,'Emit X = %10.5e\n',e_emx_sqrt);
fprintf(fileID,'Emit Y = %10.5e\n',e_emy_sqrt);
fprintf(fileID,'Electron beam normalized Emittances:\n');
fprintf(fileID,'Emit X n = %10.5e\n',norm_em_x);
fprintf(fileID,'Emit Y n = %10.5e\n',norm_em_y);
fprintf(fileID,'Other data:\n');
fprintf(fileID,'sigma_x = %10.5e\n',sigma_0);
fprintf(fileID,'sigma_y = %10.5e\n',sigma_0y);
fprintf(fileID,'gamma = %10.5e\n',gamma);
fprintf(fileID,'delta_gamma = %10.5e\n',delta_gamma);
fprintf(fileID,'electron sigma x0 = %10.5e\n',sigma_0);
fprintf(fileID,'waist W0 = %10.5e\n',W0);
fprintf(fileID,'laser_bandwidth = %10.5e\n',laser_bandwidth);
fprintf(fileID,'laserwl= %10.5e\n',laserwl)
fprintf(fileID,'a_0=4.3*(laserwl/W0)*sqrt(pulseE/sigt)= %10.5e\n',a_0);
fprintf(fileID,'----------------------------------------------------\n');
% fprintf(fileID,'NT1=((1/sqrt(12)).*gamma^2.*theta_max.^2).^2= %10.5e\n',NT1);
% fprintf(fileID,'NT2=(2*delta_gamma/gamma)^2= %10.5e\n',NT2);
% fprintf(fileID,'NT3=(2*norm_em_x^2/sigma_0^2)^2= %10.5e\n',NT3);
% fprintf(fileID,'NT3ax=(norm_em_x^2/sigma_0^2)^2= %10.5e\n',NT3ax);
% fprintf(fileID,'NT3ay=(norm_em_y^2/sigma_0y^2)^2= %10.5e\n',NT3ax);
% fprintf(fileID,'NT4=(laser_bandwidth)^2= %10.5e\n',NT4);
% fprintf(fileID,'NT5=((M2*laserwl)/(2*pi*W0))^4= %10.5e\n',NT5);
% fprintf(fileID,'NT6=((a_0^2/3)/(1+a_0^2/2))^2= %10.5e\n',NT6);
% fprintf(fileID,'BT7=sqrt(BT1*BT3)= %10.5e\n',BT7);
fclose(fileID);

end
