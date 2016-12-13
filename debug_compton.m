clear all; close all; clc;
make_path
%     stop
start_date=datestr(now);
%for new global
global ifig save_dir save_dir_start
global rflags el
global ev cs cross
[rflags] = flags_for_run;
[el] = initial_electron_beam_param;

ifig=1;

%% main parameteres
rflags.PLOTS =1;
rflags.use_ideal_beams=0;
rflags.compton=1;
%% Beams are by default  at L=0

for L=[0];
    
    %% dir for output
    %     save_dir_start=[pwd '/t3_10K_For_L_' num2str(L) '_just_one_beam/'];mkdir(save_dir_start);
    save_dir_start=[pwd '/try_Compton/'];mkdir(save_dir_start);
  
    
    %% Main loop
    E_3_full=[]; E_4_full=[]; theta_3_full=[]; theta_4_full=[]; phi_3_full=[]; phi_4_full=[];pair_full=[];
    for n_bin=[11]%[5 7 11 21 51 71 81 ]%
        close all;
        ifig=1000;
        nx_bin=n_bin;
        ny_bin=n_bin;
        nz_bin=n_bin;
        
        save_dir=[save_dir_start 'nbin_x_' num2str(nx_bin) '_y_' num2str(ny_bin) '_z_' num2str(nz_bin) '/'];mkdir(save_dir)
        

        final_plot_for_compton;
      
    end
end