% one beam several events

clear all; close all; clc;
make_path
%  stop
start_date=datestr(now);
%for new global
global home_dir DIRECTORY_FOR_CAIN BASE_DIRECTORY;
global rflags

[rflags] = flags_for_run;
rflags.PLOTS =0;
just_plots=0;

for ni=1:1:7

rflags.chargebunch = 250e-12;%Charge per electrons bunch [c]
defoc_param=1;%defocusing parameter by defoult 1 is no defocusing
home_dir=[pwd '/CAIN/'];
BASE_DIRECTORY = [pwd '/el_beams_' num2str(ni) '/'];
mkdir(BASE_DIRECTORY);
    mkdir([BASE_DIRECTORY 'initial_beam_plot/']);
    DIRECTORY_FOR_CAIN = [BASE_DIRECTORY 'cain_tmp/'];
    mkdir(DIRECTORY_FOR_CAIN);

    
    if ni==1
 beam_phasespace=dlmread(['cristina_e_beams/eli_lowen_oned_WP_ref_up_short.out.sdds.asci'],'',10,0);
    elseif ni==2
  beam_phasespace=dlmread(['cristina_e_beams/eli_lowen_oned_WP_ref_up_180MeV.out.sdds.asci'],'',10,0);
elseif ni==3
  beam_phasespace=dlmread(['cristina_e_beams/eli_lowen_oned_WP_newlayout_track_up_new_check_newsol_check.w5.asci'],'',35,0);
elseif ni==4
  beam_phasespace=dlmread(['cristina_e_beams/eli_highen_oned_WP_ref_up_short.out.sdds.asci'],'',10,0);
elseif ni==5
  beam_phasespace=dlmread(['cristina_e_beams/eli_highen_oned_WP_newlayout_track_up_new_check_newsol_rec_720MeV.out.asci'],'',10,0);
elseif ni==6
  beam_phasespace=dlmread(['cristina_e_beams/eli_highen_oned_WP_newlayout_track_up_new_check_newsol_rec_600MeV.out.asci'],'',10,0);
elseif ni==7
  beam_phasespace=dlmread(['cristina_e_beams/eli_highen_oned_WP_newlayout_track_up_new_check_newsol_rec_400MeV.out.asci'],'',10,0);
    end
[beam_property]=formating_beam_for_cain(beam_phasespace,1);
[nothing] = start_cain(beam_property,1);
rflags.PLOTS =1;
plot_expdat(BASE_DIRECTORY)

end