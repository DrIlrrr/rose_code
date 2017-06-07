function [phot_emit]=photons_emit_in_bandwidth(out_folder)

close all;
% make_path
ifig=1;

%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss

% photons_data=load(['dif__hi_600_defocus_1.5_PulseE_400mJ/photon_data_plots/photons_data.dat']);%dlmread([out_folder 'cain_tmp/cain_output_photons_' num2str(ni) '.dat'],'',1,0);%read electrons from cain
photons_data=load([out_folder 'photon_data_plots/photons_data.dat']);%dlmread([out_folder 'cain_tmp/cain_output_photons_' num2str(ni) '.dat'],'',1,0);%read electrons from cain

weigth=photons_data(1,3);
full_spectrum=photons_data(:,8)./1e3;
phot_angle=sign(photons_data(:,9)).*atan(sqrt(photons_data(:,9).^2+photons_data(:,10).^2)./photons_data(:,11));

bandwith_cm=[];
num_phot_th=[];
theta_angle=5e-7;
diapason=(0:1:101);
qq=0; num_phot_th=[]; bandwith_cm=[]; bandwith_non_norm=[];
for ni=diapason;
    qq=qq+1;
    el_angel_1=ni*theta_angle;
    aa_1=find(abs(phot_angle)<el_angel_1);    
    num_phot_th(qq)=length(full_spectrum(aa_1))*weigth;
    bandwith_cm(qq)=std(full_spectrum(aa_1))/mean(full_spectrum(aa_1));
    bandwith_non_norm(qq)=std(full_spectrum(aa_1));
    if bandwith_cm(qq)>=0.001 && bandwith_cm(qq)<=0.0052 % find angle for given bandwith
        el_angel=el_angel_1
        %         stop;
    end
end
aa=find(abs(phot_angle)<el_angel);
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
% filename = [ out_folder 'photons_plot_' num2str(ifig) ];
% fname = [ filename '.png'];
% print('-dpng', fname);


















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



phasespace(1,:)=(photons_data(aa,X_COORDINATE))';

phasespace(2,:)=(photons_data(aa,X_MOMENTUM)./(photons_data(aa,S_MOMENTUM)))';

phasespace(3,:)=(photons_data(aa,Y_COORDINATE))';

phasespace(4,:)=(photons_data(aa,Y_MOMENTUM)./(photons_data(aa,S_MOMENTUM)))';

phasespace(5,:)=photons_data(aa,TIME_COORDINATE)';

phasespace(6,:)=((photons_data(aa,ENERGY_OF_PARTICLE)-mean(photons_data(aa,ENERGY_OF_PARTICLE)))/mean(photons_data(aa,ENERGY_OF_PARTICLE)))';




% std
sig =cov(phasespace(:,:)); % Sigma matrix 6*6
b_sig =sqrt(diag(sig));     % std of the 6 variables


% emittances
i=1;emitx=sqrt(det(sig(i:i+1,i:i+1)));
i=3;emity=sqrt(det(sig(i:i+1,i:i+1)));
i=5;emits=sqrt(det(sig(i:i+1,i:i+1)));
ph_emit_in_bandwidth=[emitx ; emity ; emits];% non normalized [m rad] emittance of the 3 subspaces


fprintf('Photon Emittances in bandwidth:\n')
fprintf('Emit X = %10.5e\n',emitx)
fprintf('Emit Y = %10.5e\n',emity)
fprintf('Emit S = %10.5e\n',emits)