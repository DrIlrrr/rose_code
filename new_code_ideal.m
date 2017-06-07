clear all; close all; clc;
make_path

%   stop
start_date=datestr(now);
%for new global
global ifig save_dir save_dir_start
global rflags el
global ev cs cross gvar
[rflags] = flags_for_run;
[el] = initial_electron_beam_param;
[gvar] = global_variable;
% ifig=1;

%% main parameteres
rflags.PLOTS = 0;
rflags.use_ideal_beams=0;

%% Beams are by default  at L=0
% L=1e-2;
for scan_var=[1000 2000 3000];
    rflags.name_of_the_cain_sim=['ideal_copy_real_' num2str(scan_var) '_GeV_ES_0.1_norm_em_0.0007'];
    ifig=1;
    gvar.sigma_cut_at_end=2;
%     gvar.el_energy=250;
    %% dir for output
    save_dir_start=[pwd '/out_muons_' rflags.name_of_the_cain_sim 'sigma_cut_' num2str(gvar.sigma_cut_at_end) '/'];mkdir(save_dir_start);
    %% Load or create beams
    if rflags.use_ideal_beams==0
%         rflags.low_energy_gamma_gamma=0;
        rflags.e_gamma=1;%
        L=0
        %         [beam_1,beam_2]=FEL_gamma_gamma;
                  [beam_1,beam_2]=e_gammas_beams(L);
        %         [beam_1,beam_2]=compton(L);
        %         [beam_1,beam_2]=gamma_gamma_beams(L);
        %         [beam_1,beam_2]=e_e_beams(L);

    else
        [beam_1,beam_2]=create_ideal_beams;
    end
    %% Turn beam on angle alpha
    %     alpha=0;
    %     b2=beam_2;
    %     [beam_2]=turn_beam(alpha,b2);
    %     b2=[];
    %     beam_stat('beam_2_turned',beam_2)

    %% load crossection for gamma-gamma
    load_crossection_for_gamma_gamma;

    %% Found a limits for greed
    [x_max,x_min,y_max,y_min,z_max,z_min,l_max]=findind_min_max_of_grid(beam_1,beam_2);

    %% Plot for cheking sizes after l_max
    Plot_for_cheking_sizes_after_l_max(beam_1,beam_2,l_max)

    %% Main loop
    E_3_full=[]; E_4_full=[]; theta_3_full=[]; theta_4_full=[]; phi_3_full=[]; phi_4_full=[];pair_full=[];
    for n_bin=[51]%[5 7 11 21 51 71 81 ]%
        close all;
        ifig=1000;
        nx_bin=n_bin;
        ny_bin=n_bin;
        nz_bin=n_bin;

        save_dir=[save_dir_start 'nbin_x_' num2str(nx_bin) '_y_' num2str(ny_bin) '_z_' num2str(nz_bin) '/'];mkdir(save_dir)

        %% Make correct bining on z axes z_max-z_min is twice a l_max mesh on z are
        %  not in agriament with time step!!!!!!!!!!!!!!!!!!!
        qq=0;NumP_1=[];NumP_2=[];Ecm_tot=[];gammacm_tot=[];Ecm_step=[];pcm_tot=[];
        delta_x=[];delta_y=[];delta_z=[];


        for l=0:2*l_max/nz_bin:l_max
            qq=qq+1
            close all;
            pair_step=[];beam_1l=[];beam_2l=[];
            [beam_1l]=beam_drift(beam_1,l);
            [beam_2l]=beam_drift(beam_2,-l);
            %% main function
            [p_cm_step,gamma_cm_step,Ecm_step,NumP_1,NumP_2,delta_x,delta_y,delta_z,pair_step]=...
                running_inside_greed(beam_1l,beam_2l,x_max,x_min,nx_bin,y_max,y_min,ny_bin,z_max,z_min,nz_bin);


            pair_full=[pair_step; pair_full];
            Ecm_tot=[Ecm_tot; Ecm_step];
            gammacm_tot=[gammacm_tot; gamma_cm_step];
            pcm_tot=[pcm_tot; p_cm_step];

            visualisation_of_propagation(beam_1l,beam_2l,x_max,x_min,y_max,y_min,z_max,z_min,delta_x,delta_y,delta_z,Ecm_step,Ecm_tot,gamma_cm_step,nx_bin,ny_bin,nz_bin,qq)
            %      visualisation_of_propagation_v2(beam_1l,beam_2l,x_max,x_min,y_max,y_min,z_max,z_min,delta_x,delta_y,delta_z,Ecm_step)
        end

        %pair_info=[Vq E_3 E_4 theta_3 theta_4 phi_3 phi_4 cos_alpha cos_theta Ecm_pair gamma_cm_pair];
        save( [save_dir 'main_out_put.dat'],'pair_full','-mat')

        weight_1=beam_1(2,3);weight_2=beam_2(2,3);

        save( [save_dir 'out_put.dat'],'Ecm_tot','gammacm_tot','delta_x','delta_y','delta_z','weight_1','weight_2','-mat')


        final_plot_for_breit_wheeler
        final_plot_for_gamma_gamma;
        final_plot_for_compton;
        final_plot_for_moller;
        final_plot_for_e_gamma(gammacm_tot.*0,gammacm_tot,Ecm_tot,weight_1,weight_2,delta_x,delta_y);
        final_plot_for_TPP;
        pair_stat_plots(pair_full)

    end
end
