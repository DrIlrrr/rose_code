function [beam_1,beam_2]=gamma_gamma_beams(L)
global rflags gvar
%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss

% beam_1i_in=dlmread([pwd '/gamma_gamma_low/photons_data_5K_3sigma.dat'],'',0,0);%read electrons from cain
% beam_2i_in=dlmread([pwd '/gamma_gamma_low/photons_data_5K_3sigma.dat'],'',0,0);%read photons from cain

%  beam_1i=dlmread([pwd '/gamma_gamma_low/photons_data_10K_3sigma.dat'],'',0,0);%read electrons from cain
%  beam_2i=dlmread([pwd '/gamma_gamma_low/photons_data_10K_3sigma.dat'],'',0,0);%read photons from cain

beam_1i=dlmread([pwd '/gamma_gamma_BW_ideal/' num2str(260) 'MeVphotons_data.dat'],'',0,0);%read electrons from cain
beam_2i=dlmread([pwd '/gamma_gamma_BW_ideal/' num2str(260) 'MeVphotons_data.dat'],'',0,0);%read photons from cain

 % smaller_1_119999=find(beam_1i_in(:,8)<1.19e6);% crossection for gamma-gamma is just for 2.4MeV in CoM
 % smaller_2_119999=find(beam_2i_in(:,8)<1.18e6);
 %
 % beam_1i=beam_1i_in(smaller_1_119999,:);
 % beam_2i=beam_2i_in(smaller_2_119999,:);


[beam_1]=beam_drift(beam_1i,L);
[beam_2]=beam_drift(beam_2i,L);

beam_1(:,4)=-beam_1(:,7);
beam_2(:,4)=-beam_2(:,7);


%% create a second beam as the reflection of first one
% beam_2=beam_1;
beam_2(:,4)=-beam_2(:,4);
% make second beam fully mirrored
beam_2(:,9)=-beam_2(:,9);
beam_2(:,10)=-beam_2(:,10);
beam_2(:,5)=-beam_2(:,5);
beam_2(:,6)=-beam_2(:,6);


%% beam stat!!!!!!!!!!!!!!!!!!

beam_stat('beam_1_initial',beam_1)
beam_stat('beam_2_initial',beam_2)


%%
beam_2(:,11)=-beam_2(:,11);% thay a propagate face to face

% %% put beams around z=0
% beam_1(:,4)=beam_1(:,4)-max(beam_1(:,4));
% beam_2(:,4)=beam_2(:,4)-min(beam_2(:,4));


[x_max_temp,x_min_temp,y_max_temp,y_min_temp,z_max_temp,z_min_temp,l_max_temp]=findind_min_max_of_grid(beam_1,beam_2);

% lim_cut=1*std(beam_1(:,5));
% size(beam_1)
% size(beam_2)
% [beam_1,beam_2]=cut_tails_at_begin(beam_1,beam_2,lim_cut);
[beam_1,beam_2]=cut_tails_at_end(beam_1,beam_2,l_max_temp);
% size(beam_1)
% size(beam_2)
% stop

%% put beams around z=0
beam_1(:,4)=beam_1(:,4)-max(beam_1(:,4));
beam_2(:,4)=beam_2(:,4)-min(beam_2(:,4));

%% beam stat!!!!!!!!!!!!!!!!!!

beam_stat('beam_1_used',beam_1)
beam_stat('beam_2_used',beam_2)

% %% cut tails in x-y plaine
%
%
%
% in_10mum=find(abs(beam_1(:,5))<0.7e-5 & abs(beam_1(:,6))<0.7e-5);
% beam_temp=beam_1(in_10mum,:);
% beam_1=[];
% beam_1=beam_temp;
% beam_temp=[];in_10mum=[];
%
%
% in_10mum=find(abs(beam_2(:,5))<0.7e-5 & abs(beam_2(:,6))<0.7e-5);
% size(beam_2)
% beam_temp=beam_2(in_10mum,:);
% beam_2=[];
% beam_2=beam_temp;
%
% beam_temp=[];in_10mum=[];
