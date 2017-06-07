% one beam several events

clear all; close all; clc;
make_path
%      stop
start_date=datestr(now);
%for new global
global home_dir DIRECTORY_FOR_CAIN BASE_DIRECTORY;
global rflags el
[rflags] = flags_for_run;
[el] = initial_electron_beam_param;
home_dir=[pwd '/CAIN/'];

rflags.PLOTS =1;
just_plots=0;

for var_for_scan=[1 2 3];

    el.NUMBER_OF_MACROPARTICLES=1e5; % Number of macroparticles

    %% electron beam
    el.chargebunch = 20e-9;%Charge per electrons bunch [c] pico->10^-12
    el.initial_beam_energy_MeV=var_for_scan*1e3; % initial energy in [MeV]
    el.energy_spread_initial=1e-1;%  initial relative energy spread (not in [%])
    el.sigma_e_x=2.75*1e-6; % IP vertical electron beam size [m]
    el.sigma_e_y=2.75*1e-6; % IP horizontal electron beam size [m]
    el.norm_emit_x=7*1e-4; %Normilized emittance x [m rad]
    el.norm_emit_y=7*1e-4; %Normilized emittance y [m rad]

    %% laser
    rflags.pulseE=0.5; %laser puse energy [J]
    rflags.sigLr=5;%/2; % given in [mu m] micro meter like 2 weist w0=28;
    rflags.laserwl=800; % laser wavelenth [nm] nano meters
    rflags.sigt=5; %pulse length [ps]

    name=['test_ideal_copy_real_' num2str(el.initial_beam_energy_MeV) '_GeV_ES_' num2str(el.energy_spread_initial) '_norm_em_' num2str(el.norm_emit_x) ];
    BASE_DIRECTORY = [pwd '/' name '/'];

    %%
    mkdir(BASE_DIRECTORY);
    mkdir([BASE_DIRECTORY 'initial_beam_plot/']);
    DIRECTORY_FOR_CAIN = [BASE_DIRECTORY 'cain_tmp/'];
    mkdir(DIRECTORY_FOR_CAIN);
    %%
    %     [beam_phasespace] = defocusing_beam(beam_phasespace,defoc_param);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[beam_property]=formating_beam_for_cain(beam_phasespace,1);
    if just_plots==0
        for turn_number=1:1:10
            turn_number

            if turn_number==1
                [beam_property] = electron_beam_initial;
                %                 [beam_property] = electron_beam_initial_banana_shape;
                %                 [beam_property]=formating_beam_for_cain(beam_phasespace,1);
                [b_emit] = beam_emit(beam_property);
            end

            [nothing] = start_cain(beam_property,turn_number);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%
        end% end for turn_number
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        all_in_on_function(BASE_DIRECTORY,turn_number);
    end
    plot_expdat(BASE_DIRECTORY);
    photons_plots_simple(BASE_DIRECTORY);
    %data_for_formula(BASE_DIRECTORY,name);
    %an_formula(name);
    %     photons_plots(BASE_DIRECTORY)
    %     out_all_stat=[pwd '/WP_' num2str(var_for_scan) '/'];
    %     plot_all_in_one(BASE_DIRECTORY,out_all_stat)

end
