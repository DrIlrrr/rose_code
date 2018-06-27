
clear all; close all; clc;
make_path
ifig=1;
%%
beam_1=dlmread(['/Users/drilrrr/Desktop/nicola/beams/ideal_beam/Qe_19300Pc_pulseE_1_STOKES_0_0_1/cain_tmp/exp.dat'],'',0,0);%read electrons from cain
beam_2=dlmread(['/Users/drilrrr/Desktop/nicola/beams/ideal_beam/Qe_19300Pc_pulseE_1_STOKES_0_0_1/photon_data_plots/photons_data.dat'],'',0,0);%read photons from cain

    %% Found a limits for greed
    [x_max,x_min,y_max,y_min,z_max,z_min,l_max]=findind_min_max_of_grid(beam_1,beam_2);

%%
nz_bin=21;
nx_bin=21;
ny_bin=21;
delta_x=(x_max-x_min)/nx_bin;
delta_y=(y_max-y_min)/ny_bin;
delta_z=(z_max-z_min)/nz_bin;

for L=0:2*l_max/nz_bin:l_max
[beam_1]=beam_drift(beam_1,L);



figure(ifig)
ifig=ifig+1;
scatter3(beam_1(:,6),beam_1(:,4),beam_1(:,5),'.b')

%to have a real delta_xyz on a plot as grid
     ax = gca;
     ax.XTick = [y_min:delta_y:y_max];
     ay = gca;
     ay.YTick = [z_min:delta_z:z_max];
     az = gca;
     az.ZTick = [x_min:delta_x:x_max];
zlim([x_min x_max])
xlim([y_min y_max])
ylim([z_min z_max])
view(-90,30);
end

