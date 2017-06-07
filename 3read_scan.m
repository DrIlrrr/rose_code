clear all; close all; clc;
make_path
%     stop
start_date=datestr(now);
%for new global
global ifig save_dir save_dir_start
global rflags el
[rflags] = flags_for_run;
[el] = initial_electron_beam_param;

ifig=1;



out_dir_start=[pwd '/TO_DELETE_temp_out_e_gamma/'];mkdir(save_dir);


qq=0;
for n_bin=[21]% 101]
    qq=qq+1;
    nx_bin=n_bin;
    ny_bin=n_bin;
    nz_bin=n_bin;
    
    out_dir=[out_dir_start 'nbin_x_' num2str(nx_bin) '_y_' num2str(ny_bin) '_z_' num2str(nz_bin) '/'];
    
    load( [out_dir 'out_put.dat'],'Ecm_tot','gammacm_tot','-mat')
    
    Ecm_tot_std(qq)=std(Ecm_tot);
    
    
     figure(ifig)
    ifig=ifig+1;
    subplot (2,1,1)
    hold on
    hist(Ecm_tot,20)
    plot
    hold off
    grid on
    xlabel('E cm total')
    subplot (2,1,2)
    hist(gammacm_tot,20)
    grid on
    xlabel('\gamma cm total')
    suptitle(['bin x ' num2str(nx_bin) ' y ' num2str(ny_bin) ' z ' num2str(nz_bin) ])
    filename = [ 'gamma_EXP_fig_' num2str(ifig)];
    fname = [ filename '.png'];
    print('-dpng', fname);
    
    
end


figure(ifig)
ifig=ifig+1;
plot([5 7 11 21 51 71 81],Ecm_tot_std,'-o')
grid on
xlabel('n bin')
ylabel('std E_{cm}')
filename = [save_dir_start 'compare_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);
