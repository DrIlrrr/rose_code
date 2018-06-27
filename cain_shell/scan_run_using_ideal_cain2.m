% one beam several events

clear all; close all; clc;
make_path
%      stop
start_date=datestr(now);
%for new global
global home_dir DIRECTORY_FOR_CAIN BASE_DIRECTORY;
global rflags el
[rflags] = flags_for_run;
[el] = initial_electron_beam_param;
home_dir=[pwd '/CAIN/'];

rflags.PLOTS =1;
just_plots=0;
qq=0;
sigma_x_var=[1 1.1 1.2 1.3 1.4 1.5 2 3 4];

for var_for_scan =sigma_x_var;
qq=qq+1;

    el.NUMBER_OF_MACROPARTICLES=1e4; % Number of macroparticles
    %% electron beam
    el.chargebunch = 22e-9;%Charge per electrons bunch [c] pico->10^-12
    el.initial_beam_energy_MeV=1600; % initial energy in [MeV]
    el.energy_spread_initial=0.1;%  initial relative energy spread (not in [%])
    el.sigma_e_x=var_for_scan*1e-6; % IP vertical electron beam size [m]
    el.sigma_e_y=var_for_scan*1e-6; % IP horizontal electron beam size [m]
    el.norm_emit_x=1e-5; %Normilized emittance x [m rad]
    el.norm_emit_y=1e-5; %Normilized emittance y [m rad]
    el.bunch_length_initial=5e-6; % intial bunch length [m]
    %% laser
    rflags.pulseE=5; %laser pulse energy [J]
    rflags.sigLr=5;%/2; % given in [mu m] micro meter like 2 weist w0=28;
    rflags.laserwl=800; % laser wavelenth [nm] nano meters
    rflags.sigt=3; %pulse length [ps]

    name=['ideal_sx_' num2str(el.sigma_e_x*1e6) '_mum_' num2str(el.initial_beam_energy_MeV/1e3) '_GeV_ES_' num2str(el.energy_spread_initial) '_norm_em_' num2str(el.norm_emit_x) ];
    BASE_DIRECTORY = [pwd '/' name '/'];

beam_in=[];
beam_in=load([BASE_DIRECTORY 'cain_tmp/exp.dat']);
%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss
% WEIGHT=3;
TIME_COORDINATE=4;%now we start use s(m) for caine
X_COORDINATE=5;
Y_COORDINATE=6;
ENERGY_OF_PARTICLE=8;
X_MOMENTUM=9;
Y_MOMENTUM=10;
S_MOMENTUM=11;

PX(qq)=std(1e-6*beam_in(:,X_MOMENTUM));
PY(qq)=std(1e-6*beam_in(:,Y_MOMENTUM));
%%
e_x=(beam_in(:,X_COORDINATE));
e_xp=(beam_in(:,X_MOMENTUM)./(beam_in(:,S_MOMENTUM)));
e_y=(beam_in(:,Y_COORDINATE));
e_yp=(beam_in(:,Y_MOMENTUM)./(beam_in(:,S_MOMENTUM)));
e_s=beam_in(:,TIME_COORDINATE)';
e_ens=((beam_in(:,ENERGY_OF_PARTICLE)-mean(beam_in(:,ENERGY_OF_PARTICLE)))/mean(beam_in(:,ENERGY_OF_PARTICLE)))';
e_E=beam_in(:,ENERGY_OF_PARTICLE);

mean_e_Px=mean(beam_in(:,X_MOMENTUM));
mean_e_Py=mean(beam_in(:,Y_MOMENTUM));
std_e_Px=std(beam_in(:,X_MOMENTUM));
std_e_Py=std(beam_in(:,Y_MOMENTUM));

e_emx_sqrt(qq)=sqrt(mean((e_x-mean(e_x)).^2)*mean((e_xp-mean(e_xp)).^2)-mean((e_x-mean(e_x)).*(e_xp-mean(e_xp)))^2);
e_emy_sqrt(qq)=sqrt(mean((e_y-mean(e_y)).^2)*mean((e_yp-mean(e_yp)).^2)-mean((e_y-mean(e_y)).*(e_yp-mean(e_yp)))^2);




photons_data=[];
photons_data=load([BASE_DIRECTORY 'photon_data_plots/photons_data.dat']);
  
weigth=photons_data(1,3);
full_spectrum=photons_data(:,8)./1e3;
number_of_photons(qq)=length(full_spectrum)*weigth;  

end

ifig=1

figure(ifig)
ifig=ifig+1;
hold on
plot(sigma_x_var,e_emx_sqrt,'-ob')
plot(sigma_x_var,e_emy_sqrt,'-sb')
hold off
grid on
set(gca,'FontSize',16)
legend('\varepsilon_x','\varepsilon_y',0)
ylabel('emittance')
xlabel('\sigma_x')
filename = [ 'c_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);



figure(ifig)
ifig=ifig+1;
hold on
plot(sigma_x_var,PX,'-ob')
plot(sigma_x_var,PY,'-sb')
hold off
grid on
set(gca,'FontSize',16)
legend('P_x','P_y',0)
ylabel('transvers momentum')
xlabel('\sigma_x')
filename = [ 'c_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);



figure(ifig)
ifig=ifig+1;
plot(sigma_x_var,number_of_photons,'-ob')
grid on
set(gca,'FontSize',16)
legend('P_x','P_y',0)
ylabel('# Compton photons')
xlabel('\sigma_x')
filename = [ 'c_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);