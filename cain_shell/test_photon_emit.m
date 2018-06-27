clear all; close all; clc;
make_path


%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss

photons_data=load(['dif__hi_600_defocus_1.5_PulseE_400mJ/photon_data_plots/photons_data.dat']);%dlmread([out_folder 'cain_tmp/cain_output_photons_' num2str(ni) '.dat'],'',1,0);%read electrons from cain


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



phasespace(1,:)=(photons_data(:,X_COORDINATE))';

phasespace(2,:)=(photons_data(:,X_MOMENTUM)./(photons_data(:,S_MOMENTUM)))';

phasespace(3,:)=(photons_data(:,Y_COORDINATE))';

phasespace(4,:)=(photons_data(:,Y_MOMENTUM)./(photons_data(:,S_MOMENTUM)))';

phasespace(5,:)=photons_data(:,TIME_COORDINATE)';

phasespace(6,:)=((photons_data(:,ENERGY_OF_PARTICLE)-mean(photons_data(:,ENERGY_OF_PARTICLE)))/mean(photons_data(:,ENERGY_OF_PARTICLE)))';




% std
sig =cov(phasespace(:,:)); % Sigma matrix 6*6
b_sig =sqrt(diag(sig));     % std of the 6 variables


% emittances
i=1;emitx=sqrt(det(sig(i:i+1,i:i+1)));
i=3;emity=sqrt(det(sig(i:i+1,i:i+1)));
i=5;emits=sqrt(det(sig(i:i+1,i:i+1)));
phot_emit=[emitx ; emity ; emits];% non normalized [m rad] emittance of the 3 subspaces 


fprintf('Photon Emittances:\n')
fprintf('Emit X = %10.5e\n',emitx)
fprintf('Emit Y = %10.5e\n',emity)
fprintf('Emit S = %10.5e\n',emits)



