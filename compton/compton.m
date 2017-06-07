function [beam_1,beam_2]=compton(L)
global rflags el
%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss

rflags.compton=1;

%% create a electron beam
% % el.NUMBER_OF_MACROPARTICLES=1e3;
% [beam_1] = electron_beam_initial;
% 
%% create a laser beam
% [beam_2] = create_laser_beam;

%% load beam
beam_1i=dlmread([pwd '/gamma_gamma_BW_ideal/' num2str(260) 'MeVelectrons_data.dat'],'',1,0);%read electrons from cain
beam_2i=dlmread([pwd '/gamma_gamma_BW_ideal/' num2str(260) 'MeVphotons_data.dat'],'',0,0);%read photons from cain

[beam_1]=beam_drift(beam_1i,L);
[beam_2]=beam_drift(beam_2i,L);
%%
 beam_1(:,4)=-beam_1(:,7);
 beam_2(:,4)=-beam_2(:,7);

%% create a second beam as the reflection of first one
beam_2(:,4)=-beam_2(:,4);
% make second beam fully mirrored
beam_2(:,9)=-beam_2(:,9);
beam_2(:,10)=-beam_2(:,10);
beam_2(:,5)=-beam_2(:,5);
beam_2(:,6)=-beam_2(:,6);

%% put beams around z=0
beam_1(:,4)=beam_1(:,4)-max(beam_1(:,4));
beam_2(:,4)=beam_2(:,4)-min(beam_2(:,4));

%%
beam_2(:,11)=-beam_2(:,11);% thay a propagate face to face

%% beam stat!!!!!!!!!!!!!!!!!!

beam_stat('electrons_used',beam_1)
beam_stat('laser_used',beam_2)
