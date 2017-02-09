% one beam several events

clear all; close all; clc;
make_path
     %  stop
start_date=datestr(now);
%for new global
global home_dir DIRECTORY_FOR_CAIN BASE_DIRECTORY;
global rflags
[rflags] = flags_for_run;
home_dir=[pwd '/CAIN/'];

rflags.PLOTS =1;
just_plots=0;


for var_for_scan=1:1:1;


    rflags.chargebunch = 19.3e-9;%Charge per electrons bunch [c]
    defoc_param=1;%defocusing parameter by defoult 1 is no defocusing

    %     name_n=['Qe_' num2str(rflags.chargebunch*1e12) 'Pc_pulseE_' num2str(rflags.pulseE) '_energy_factor_' num2str(var_for_scan)];
    %     BASE_DIRECTORY = [pwd '/' name_n '/'];
    beam_phasespace=dlmread(['nicola_beam_24_01_2016/005/Phase_space_09.dat'],'',0,0);

    beam_phasespace(:,5)=beam_phasespace(:,5)./3e8;
    beam_phasespace(:,6)=6.*beam_phasespace(:,6)./1e6;


    name_n=['Qe_1_Pc_pulseE_' num2str(rflags.pulseE) '_energy_*6_500MeV_cut'];
    BASE_DIRECTORY = [pwd '/' name_n '/'];


    mkdir(BASE_DIRECTORY);
    mkdir([BASE_DIRECTORY 'initial_beam_plot/']);
    DIRECTORY_FOR_CAIN = [BASE_DIRECTORY 'cain_tmp/'];
    mkdir(DIRECTORY_FOR_CAIN);

    %     beam_phasespace(:,6)=beam_phasespace(:,6)-114.7/0.511;
    %     defoc_param=1.7;
    %     [beam_phasespace] = defocusing_beam(beam_phasespace,defoc_param);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [beam_property]=formating_beam_for_cain(beam_phasespace,1);


    cut_beam=find(beam_property(8,:)>500e6);

    newbeam=beam_property(:,cut_beam);

    beam_property=[];
    beam_property=newbeam;
    newbeam=[];


  %  figure(1)
  %  hist(beam_property(8,:),20);

    %     stop

    if just_plots==0
        for turn_number=1:1:1
            turn_number

            %             if turn_number==1
            %                 [beam_property]=formating_beam_for_cain(beam_phasespace,1);
            %                 %                 stop
            %                 [b_emit] = beam_emit(beam_property);
            %             end

            [nothing] = start_cain(beam_property,turn_number);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%
        end% end for turn_number
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        all_in_on_function(BASE_DIRECTORY,turn_number);
    end
    plot_expdat(BASE_DIRECTORY)
    photons_plots(BASE_DIRECTORY)
    for ni=1:1:turn_number
        system(['rm ' BASE_DIRECTORY 'cain_tmp/cain_output_photons_' num2str(ni) '.dat']); % delete exp.dat
    end
end
