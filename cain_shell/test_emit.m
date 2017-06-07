% one beam several events

clear all; close all; clc;
make_path

start_date=datestr(now);
%for new global
global home_dir DIRECTORY_FOR_CAIN BASE_DIRECTORY;
global rflags
[rflags] = flags_for_run;
home_dir=[pwd '/CAIN/'];

rflags.PLOTS =1;
just_plots=0;


   
    
    
%     rflags.STOKES=[0 1 0];%STOKES(var_for_scan,:);
%     
%     
%     rflags.pulseE=1; %laser puse energy [J]
%     rflags.angle=0; %initial scattered angle [rad]
%     rflags.laserwl=1000; % laser wavelenth [nm] nano meters
%     rflags.sigLr=5; % given in [mu m] micro meter like 2 weist w0=28;
%     rflags.sigt=3; %pulse length [ps]
    
    
    defoc_param=1;%defocusing parameter by defoult 1 is no defocusing
    BASE_DIRECTORY = [pwd '/11dif_PulseE_' num2str(rflags.pulseE) 'J_defocus_' num2str(defoc_param) '/'];
    
    
%       beam_phasespace=dlmread(['Alberto_beam/Tstep_file.txt'],'',0,0);
%     beam_phasespace=dlmread(['Alberto_beam/astra00.sdds'],'',11,0);
    beam_phasespace=dlmread(['cristina_e_beams/eli_highen_oned_WP_newlayout_track_up_new_check_newsol_rec_600MeV.out.asci'],'',10,0);
      
    
    mkdir(BASE_DIRECTORY);
    mkdir([BASE_DIRECTORY 'initial_beam_plot/']);
    DIRECTORY_FOR_CAIN = [BASE_DIRECTORY 'cain_tmp/'];
    mkdir(DIRECTORY_FOR_CAIN);
    
    %  beam_phasespace(:,6)=beam_phasespace(:,6);%-114.7/0.511;
    
      [beam_phasespace] = defocusing_beam(beam_phasespace,defoc_param);
    



[beam_property]=formating_beam_for_cain(beam_phasespace,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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



phasespace(1,:)=(beam_property(X_COORDINATE,:))';

phasespace(2,:)=(beam_property(X_MOMENTUM,:)./(beam_property(S_MOMENTUM,:)))';

phasespace(3,:)=(beam_property(Y_COORDINATE,:))';

phasespace(4,:)=(beam_property(Y_MOMENTUM,:)./(beam_property(S_MOMENTUM,:)))';

phasespace(5,:)=beam_property(TIME_COORDINATE,:)';

phasespace(6,:)=((beam_property(ENERGY_OF_PARTICLE,:)-mean(beam_property(ENERGY_OF_PARTICLE,:)))/mean(beam_property(ENERGY_OF_PARTICLE,:)))';




% std
sig =cov(phasespace(:,:)'); % Sigma matrix 6*6
b_sig =sqrt(diag(sig));     % std of the 6 variables


% emittances
i=1;emitx=sqrt(det(sig(i:i+1,i:i+1)));
i=3;emity=sqrt(det(sig(i:i+1,i:i+1)));
i=5;emits=sqrt(det(sig(i:i+1,i:i+1)));
b_emit=[emitx ; emity ; emits]% non normalized [m rad] emittance of the 3 subspaces 