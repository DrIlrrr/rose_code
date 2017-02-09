function [beam_1,beam_2]=e_gammas_beams(L);
% load electron and gamma beam
global ifig rflags save_dir_start



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

% rflags.e_gamma=1;

%% load beam profiles

beam_1i=dlmread([pwd '/cain_shell/' rflags.name_of_the_cain_sim '/cain_tmp/exp.dat'],'',0,0);%read electrons from cain
beam_2i=dlmread([pwd '/cain_shell/' rflags.name_of_the_cain_sim '/photon_data_plots/photons_data.dat'],'',0,0);%read photons from cain

% temp_beam_ind=find(beam_2i(:,7)>0.02e-6);
% beam_t=beam_2i(temp_beam_ind,:);
% beam_2i=[];
% beam_2i=beam_t;


%  beam_1i(:,4)=-beam_1i(:,7);
 beam_2i(:,4)=-beam_2i(:,7);
% propagation of the beam before collision
[beam_1]=beam_drift(beam_1i,L);
[beam_2]=beam_drift(beam_2i,L);


%% create a second beam as the reflection of first one
% beam_2=beam_1;
beam_2(:,4)=-beam_2(:,4);
% make second beam fully mirrored
beam_2(:,9)=-beam_2(:,9);
beam_2(:,10)=-beam_2(:,10);
beam_2(:,5)=-beam_2(:,5);
beam_2(:,6)=-beam_2(:,6);
beam_2(:,11)=-beam_2(:,11);% thay a propagate face to face

%% beam stat!!!!!!!!!!!!!!!!!!

beam_stat('beam_1_initial',beam_1)
beam_stat('beam_2_initial',beam_2)

%% put beams around z=0
beam_1(:,4)=beam_1(:,4)-max(beam_1(:,4));
beam_2(:,4)=beam_2(:,4)-min(beam_2(:,4));


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



%
% % %% make a small cut of tails
% % beam_2=beam_2(find(abs(beam_2(:,4))<0.2e-4),:);
% % % beam_1(:,4)=beam_1(:,7);
% % beam_1=beam_1(find(abs(beam_1(:,4))<0.2e-4),:);
%
% %%
% beam_2(:,11)=-beam_2(:,11);% thay a propagate face to face
%
% %% put beams around z=0
% beam_1(:,4)=beam_1(:,4)-max(beam_1(:,4));
% beam_2(:,4)=beam_2(:,4)-min(beam_2(:,4));
%
%
% %% cut tails in x-y plaine
%
% in_10mum_start=find(abs(beam_1(:,5))<2e-5 & abs(beam_1(:,6))<2e-5);
% beam_temp=beam_1(in_10mum_start,:);
% beam_1=[];
% beam_1=beam_temp;
% beam_temp=[];in_10mum=[];
%
%
% in_10mum_start=find(abs(beam_2(:,5))<2e-5 & abs(beam_2(:,6))<2e-5);
% size(beam_2)
% beam_temp=beam_2(in_10mum_start,:);
% beam_2=[];
% beam_2=beam_temp;
%
% beam_temp=[];in_10mum=[];
%
% %% cut tails in x-y plaine after propagation on l_max
%
%
% % [x_max_temp,x_min_temp,y_max_temp,y_min_temp,z_max_temp,z_min_temp,l_max_temp]=findind_min_max_of_grid(beam_1,beam_2);
% %
% % lim_cut=0.9e-5;
% % % size(beam_1)
% % % size(beam_2)
% % [beam_1,beam_2]=cut_tails_at_end(beam_1,beam_2,0,lim_cut);
% % % size(beam_1)
% % % size(beam_2)
% % % stop
%
%
%
%
% beam_stat('beam_1_cuted_start',beam_1)
%
% beam_stat('beam_2_cuted_start',beam_2)
%
%
% %% cut tails in x-y plaine after propagation on l_max
%
%
% [x_max_temp,x_min_temp,y_max_temp,y_min_temp,z_max_temp,z_min_temp,l_max_temp]=findind_min_max_of_grid(beam_1,beam_2);
%
% lim_cut=2e-6;
% size(beam_1)
% size(beam_2)
% [beam_1,beam_2]=cut_tails_at_end(beam_1,beam_2,l_max_temp);%,lim_cut);
% size(beam_1)
% size(beam_2)
% % stop
% beam_stat('beam_1_cuted_end',beam_1)
% beam_stat('beam_2_cuted_end',beam_2)
