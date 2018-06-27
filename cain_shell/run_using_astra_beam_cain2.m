% one beam several events

clear all; close all; clc;
make_path
%     stop
start_date=datestr(now);
%for new global
global home_dir DIRECTORY_FOR_CAIN BASE_DIRECTORY;
global rflags
[rflags] = flags_for_run;
home_dir=[pwd '/CAIN/'];

rflags.PLOTS =1;
just_plots=0;


%
% STOKES(1,:)=[1 0 0];
% STOKES(2,:)=[0 1 0];
% STOKES(3,:)=[0 0 1];


for var_for_scan=[3]%1:1:1;
    
    %     rflags.STOKES=STOKES(var_for_scan,:);
    rflags.angle=0*(pi/180);
    
    rflags.pulseE=1; %laser puse energy [J]
    rflags.angle=0; %initial scattered angle [rad]
    rflags.laserwl=1000; % laser wavelenth [nm] nano meters
    rflags.sigLr=5; % given in [mu m] micro meter like 2 weist w0=28;
    rflags.sigt=3; %pulse length [ps]
    
    rflags.chargebunch = 250e-12;%Charge per electrons bunch [c]
    defoc_param=1;%defocusing parameter by defoult 1 is no defocusing
    BASE_DIRECTORY = [pwd '/10K_no_halo_sigma_' num2str(var_for_scan) '_defocus_' num2str(defoc_param) '/'];
    
    
    astra_beam=load('astra_beams/250MeV_10k.0530.002');
    beam_phasespace=convert_astra_to_elegant(astra_beam);
    
    
    %     beam_phasespace=dlmread(['marcello_beam_30_03_2016/BeamFocused.ascii'],'',11,0);
    %     beam_phasespace=dlmread(['marcello_beam_30_03_2016/BeamFocused.ascii'],'',0,0);
    
    
    
    
    mkdir(BASE_DIRECTORY);
    mkdir([BASE_DIRECTORY 'initial_beam_plot/']);
    DIRECTORY_FOR_CAIN = [BASE_DIRECTORY 'cain_tmp/'];
    mkdir(DIRECTORY_FOR_CAIN);
    
    %     beam_phasespace(:,6)=beam_phasespace(:,6)-70/0.511;
    
    [beam_phasespace] = defocusing_beam(beam_phasespace,defoc_param);
    beam_phasespace(:,1)=beam_phasespace(:,1)-mean(beam_phasespace(:,1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[beam_property]=formating_beam_for_cain(beam_phasespace,1);
    if just_plots==0
        for turn_number=1:1:10
            turn_number
            
            if turn_number==1
                [beam_property]=formating_beam_for_cain(beam_phasespace,1);
                
                %% Cut halo
                lim_cut=var_for_scan*std(beam_property(5,:));
                in_10mum_1=find((beam_property(5,:)).^2+(beam_property(6,:)).^2<lim_cut^2);
                beam_temp=beam_property(:,in_10mum_1);
                beam_property=[];
                beam_property=beam_temp;
                beam_temp=[];%in_10mum_1=[];
                
                [b_emit] = beam_emit(beam_property);
            end
            
            [nothing] = start_cain(beam_property,turn_number);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
        end% end for turn_number
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        all_in_on_function(BASE_DIRECTORY,turn_number);
    end
    plot_expdat(BASE_DIRECTORY)
    photons_plots(BASE_DIRECTORY)
    
    
    %      photons_emit_in_bandwidth(BASE_DIRECTORY)
    
    %     analitical_f(BASE_DIRECTORY,beam_property,b_emit)
    %     photons_emit(BASE_DIRECTORY)
    %     stop
    
    
    
    
end
