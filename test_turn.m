clear all; close all; clc;
make_path
%           stop
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
rflags.PLOTS =1;

ifig=1;
gvar.sigma_cut_at_end=2;

%% dir for output
save_dir_start=[pwd '/test_turn/'];mkdir(save_dir_start);
% [beam_1,beam_2]=gg_turn;
beam_1i=dlmread([pwd '/gamma_gamma_BW_ideal/260MeVphotons_data.dat'],'',0,0);%read electrons from cain
beam_2i=dlmread([pwd '/gamma_gamma_BW_ideal/260MeVphotons_data.dat'],'',0,0);%read photons from cain

L=0;
gvar.sigma_cut_at_end=10;
[beam_1,beam_2]=cut_tails_at_end(beam_1i,beam_2i,0);

[beam_1]=beam_drift(beam_1i,L);
[beam_2]=beam_drift(beam_2i,L);

beam_1(:,4)=-beam_1(:,7);
beam_2(:,4)=-beam_2(:,7);
% beam_2(:,4)=2.*beam_2(:,4);
%% create a second beam as the reflection of first one
% beam_2=beam_1;
beam_2(:,4)=-beam_2(:,4);
% make second beam fully mirrored
beam_2(:,9)=-beam_2(:,9);
beam_2(:,10)=-beam_2(:,10);
beam_2(:,5)=-beam_2(:,5);
beam_2(:,6)=-beam_2(:,6);

% %% beam stat!!!!!!!!!!!!!!!!!!
% beam_stat('beam_1_initial',beam_1)
% beam_stat('beam_2_initial',beam_2)

%%
beam_2(:,11)=-beam_2(:,11);% thay a propagate face to face




%% put beams around z=0
beam_1(:,4)=beam_1(:,4)-max(beam_1(:,4));
beam_2(:,4)=beam_2(:,4)-min(beam_2(:,4));


%%

for alpha=[0:0.01:pi];
 beam_2t=[];
  [beam_2t]=turn_beam(alpha,beam_2);

 figure(ifig)
ifig=ifig+1;
set(gca,'FontSize',16)
hold on
scatter3(beam_2(:,5),beam_2(:,6),beam_2(:,4),'.g')
scatter3(beam_2t(:,5),beam_2t(:,6),beam_2t(:,4),'.b')
hold off
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
title(['alpha=' num2str(alpha) ' [rad]; alpha*180/pi=' num2str(alpha*180/pi) ' [deg];'])
view(-90,0)
fname2 = [ 'turn_beam_' num2str(ifig) '.png'];
print('-r300','-dpng', fname2); 


end
   stop
 
%  beam_2=[];
%  beam_2=beam_2t;
 
 
figure(ifig)
ifig=ifig+1;
set(gca,'FontSize',16)
hold on
scatter3(beam_1(:,5),beam_1(:,6),beam_1(:,4),'.b')
scatter3(beam_2(:,5),beam_2(:,6),beam_2(:,4),'.r')
hold off
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')





%%
[x_max_temp,x_min_temp,y_max_temp,y_min_temp,z_max_temp,z_min_temp,l_max_temp]=findind_min_max_of_grid(beam_1,beam_2);
% [beam_1,beam_2]=cut_tails_at_end(beam_1,beam_2,l_max_temp);

n_bin=[21]
nx_bin=n_bin;
ny_bin=n_bin;
nz_bin=n_bin;

delta_x=(x_max_temp-x_min_temp)/nx_bin;
delta_y=(y_max_temp-y_min_temp)/ny_bin;
delta_z=(z_max_temp-z_min_temp)/nz_bin;
qq=0;
for l=0:2*l_max_temp/nz_bin:l_max_temp
    qq=qq+1
    
    [beam_1l]=beam_drift(beam_1,l);
    [beam_2l]=beam_drift(beam_2,-l);
    visualisation_of_propagation_v2(beam_1l,beam_2l,x_max_temp,x_min_temp,y_max_temp,y_min_temp,z_max_temp,z_min_temp,delta_x,delta_y,delta_z)
   
end
% %% beam stat!!!!!!!!!!!!!!!!!!
% beam_stat('beam_1_used',beam_1)
% beam_stat('beam_2_used',beam_2)
