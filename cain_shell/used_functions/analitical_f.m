function analitical_f(out_folder,beam_property,b_emit)
global rflags
% clear all; close all; clc;
% make_path
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
photons_data=load([out_folder 'photon_data_plots/photons_data.dat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weigth=photons_data(1,3);
full_spectrum=photons_data(:,8)./1e3;
phot_angle=sign(photons_data(:,9)).*atan(sqrt(photons_data(:,9).^2+photons_data(:,10).^2)./photons_data(:,11));

bandwith_cm=[];
num_phot_th=[];
theta_angle=5e-7;
diapason=(0:1:101);
qq=0; num_phot_th=[]; bandwith_cm=[]; 
for ni=diapason;
    qq=qq+1;
    el_angel_1=ni*theta_angle;
    aa_1=find(abs(phot_angle)<el_angel_1);    
    num_phot_th(qq)=length(full_spectrum(aa_1))*weigth;
    bandwith_cm(qq)=std(full_spectrum(aa_1))/mean(full_spectrum(aa_1));
    if bandwith_cm(qq)>=0.001 && bandwith_cm(qq)<=0.0052 % find angle for given bandwith
        el_angel=el_angel_1;
        %         stop;
    end
end
aa=find(abs(phot_angle)<el_angel);
el_angel
 bandwith=std(full_spectrum(aa))/mean(full_spectrum(aa))
% turn_number=20000;
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

phasespace(1,:)=(photons_data(aa,X_COORDINATE))';
phasespace(2,:)=(photons_data(aa,X_MOMENTUM)./(photons_data(aa,S_MOMENTUM)))';
phasespace(3,:)=(photons_data(aa,Y_COORDINATE))';
phasespace(4,:)=(photons_data(aa,Y_MOMENTUM)./(photons_data(aa,S_MOMENTUM)))';
phasespace(5,:)=photons_data(aa,TIME_COORDINATE)';
phasespace(6,:)=((photons_data(aa,ENERGY_OF_PARTICLE)-mean(photons_data(aa,ENERGY_OF_PARTICLE)))/mean(photons_data(aa,ENERGY_OF_PARTICLE)))';
% std
sig =cov(phasespace(:,:)'); % Sigma matrix 6*6
b_sig =sqrt(diag(sig));     % std of the 6 variables

% emittances
i=1;emitx=sqrt(det(sig(i:i+1,i:i+1)));
i=3;emity=sqrt(det(sig(i:i+1,i:i+1)));
i=5;emits=sqrt(det(sig(i:i+1,i:i+1)));
ph_emit_in_bandwidth=[emitx ; emity ; emits];% non normalized [m rad] emittance of the 3 subspaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_0=std(beam_property(X_COORDINATE,:));
W0=rflags.sigLr*2*1e-6;
gamma=mean(beam_property(ENERGY_OF_PARTICLE,:))/(0.511e6);
delta_gamma=std(beam_property(ENERGY_OF_PARTICLE,:)./(0.511e6));
% bandwith=0.005; look up
norm_em_x=sqrt(gamma^2-1)*b_emit(1);
laser_bandwidth=1e-3;
M2=1.2;
laserwl=rflags.laserwl*1e-9;
a_0=4.3*(laserwl/W0)*sqrt(rflags.pulseE/rflags.sigt);


first_term=(sigma_0*W0)/(gamma*sqrt(4*sigma_0^2+W0^2))

AT1=(bandwith)^2
AT2=(2*delta_gamma/gamma)^2
AT3=(2*norm_em_x^2/sigma_0^2)^2 %check rad??
AT4=(laser_bandwidth)^2
AT5=((M2*laserwl)/(2*pi*W0))^4
AT6=((a_0^2/3)/(1+a_0^2/2))^2

em_gamma_0=first_term*...
    sqrt(AT1-AT2-AT3-AT3-AT4-AT5-AT6)

emitx

emity

emits


sqrt(std(phasespace(1,:).^2)*std(phasespace(2,:).^2)-mean(phasespace(1,:).*phasespace(2,:))^2)


mean(phasespace(1,:).*phasespace(2,:))^2





















