function [beam_1,beam_2]=FEL_gamma_gamma
global rflags el gvar

%% create e electron beam
[beam_1] = FEL_gamma_initial;
% beam_stat(out_folder,beam_1); %save a beam stat

%% load beam
% beam_2=dlmread([pwd '/gamma_gamma_BW_ideal/' num2str(gvar.el_energy) 'MeVphotons_data.dat'],'',0,0);%read photons from cain
beam_2=dlmread([pwd '/gamma_gamma_low/photons_data_10K_3sigma.dat'],'',0,0);%read photons from cain


%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss
beam_2(:,8)=ones(1,length(beam_2(:,8)))*1.2e6;
beam_2(:,11)=ones(1,length(beam_2(:,11)))*1.2e6;


beam_2(:,4)=-beam_2(:,7);
%% create a second beam as the reflection of first one
% beam_2=beam_1;
 beam_2(:,11)=-beam_2(:,11);
 beam_2(:,4)=-beam_2(:,4);
% beam_2(:,5)=beam_1(:,6);
% beam_2(:,6)=beam_1(:,5);
% [beam_1l]=beam_drift(beam_1,0);
%



%% put beams around z=0
beam_1(:,4)=beam_1(:,4)-max(beam_1(:,4));
beam_2(:,4)=beam_2(:,4)-min(beam_2(:,4));

% beam_1(:,4)=beam_1(:,4).*5;
% beam_2(:,4)=beam_2(:,4)./5;

%% beam stat!!!!!!!!!!!!!!!!!!

beam_stat('beam_1_initial',beam_1)
beam_stat('beam_2_initial',beam_2)