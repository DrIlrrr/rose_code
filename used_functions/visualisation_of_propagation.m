function visualisation_of_propagation(beam_1l,beam_2l,x_max,x_min,y_max,y_min,z_max,z_min,delta_x,delta_y,delta_z,Ecm_step,Ecm_tot,gamma_cm_step,nx_bin,ny_bin,nz_bin,qq)
global ifig save_dir rflags


weight_1=beam_1l(2,3);weight_2=beam_2l(2,3);

ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end


set(gca,'FontSize',16)
subplot (2,2,[1 3])
set(gca,'FontSize',16)
hold on
%    scatter3(y_min:delta_y:y_max,z_min:delta_z:z_max,x_min:delta_x:x_max,'k')
  scatter3(beam_1l(:,6),beam_1l(:,4),beam_1l(:,5),'.b')
  scatter3(beam_2l(:,6),beam_2l(:,4),beam_2l(:,5),'.r')
hold off
grid on
%to have a real delta_xyz on a plot as grid
%     ax = gca;
%     ax.XTick = [y_min:delta_y:y_max];
%     ay = gca;
%     ay.YTick = [z_min:delta_z:z_max];
%     az = gca;
%     az.ZTick = [x_min:delta_x:x_max];
zlim([x_min x_max])
xlim([y_min y_max])
ylim([z_min z_max])
view(-55,15);
xlabel('Y')
ylabel('Z')
zlabel('X')
 if rflags.TPP==1
     title('blu-e_{-}, red-\gamma')
% legend('e_{-}','\gamma','Location','northwest')
 end
subplot (2,2,2)
set(gca,'FontSize',16)
% nbin_plot=20;%floor(max(full_spectrum(aa))-min(full_spectrum(aa)));% nbin to make spectrum per KeV
% xs_step=linspace(min(Ecm_step),max(Ecm_step),nbin_plot);
% bar(xs_step,hist(Ecm_step,nbin_plot)*weight_1*weight_2,'grouped')%,'hist','g')
hist(Ecm_step,20)
grid on
xlabel('E cm step')
subplot (2,2,4)
set(gca,'FontSize',16)
% xs_tot=linspace(min(Ecm_tot),max(Ecm_tot),nbin_plot);
% bar(xs_tot,hist(Ecm_tot,nbin_plot)*weight_1*weight_2,'grouped')%,'hist','g')
hist(Ecm_tot,20)
grid on
xlabel('E cm total')
suptitle(['bin x ' num2str(nx_bin) ' y ' num2str(ny_bin) ' z ' num2str(nz_bin) ' step ' num2str(qq)])
filename = [save_dir 'EXP_fig_' num2str(ifig)];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
% print('-depsc', fname);
print('-r300','-dpng', fname2);


%     figure(ifig)
%     ifig=ifig+1;
%     subplot (2,2,[1 3])
%     hold on
%     %    scatter3(y_min:delta_y:y_max,z_min:delta_z:z_max,x_min:delta_x:x_max,'k')
%     scatter3(beam_1l(:,6),beam_1l(:,4),beam_1l(:,5),'.b')
%     scatter3(beam_2l(:,6),beam_2l(:,4),beam_2l(:,5),'.r')
%     hold off
%     grid on
% %to have a real delta_xyz on a plot as grid
%     %     ax = gca;
% %     ax.XTick = [y_min:delta_y:y_max];
% %     ay = gca;
% %     ay.YTick = [z_min:delta_z:z_max];
% %     az = gca;
% %     az.ZTick = [x_min:delta_x:x_max];
%     zlim([x_min x_max])
%     xlim([y_min y_max])
%     ylim([z_min z_max])
%     view(-55,15);
%     xlabel('Y')
%     ylabel('Z')
%     zlabel('X')
%     subplot (2,2,2)
%     hist(Ecm_step)
%     grid on
%     xlabel('E cm step')
%      subplot (2,2,4)
%     hist(Ecm_tot)
%     grid on
%     xlabel('E cm total')
%     suptitle(['bin x ' num2str(nx_bin) ' y ' num2str(ny_bin) ' z ' num2str(nz_bin) ' step ' num2str(qq)])
%     filename = [save_dir 'EXP_fig_' num2str(ifig)];
%     fname = [ filename '.png'];
%     print('-dpng', fname);
%
%
%
%      figure(ifig)
%     ifig=ifig+1;
%     subplot (2,3,[1 4])
%     hold on
%     %    scatter3(y_min:delta_y:y_max,z_min:delta_z:z_max,x_min:delta_x:x_max,'k')
%     scatter3(beam_1l(:,6),beam_1l(:,4),beam_1l(:,5),'.b')
%     scatter3(beam_2l(:,6),beam_2l(:,4),beam_2l(:,5),'.r')
%     hold off
%     grid on
%     ax = gca;
%     ax.XTick = [y_min:delta_y:y_max];
%     ay = gca;
%     ay.YTick = [z_min:delta_z:z_max];
%     az = gca;
%     az.ZTick = [x_min:delta_x:x_max];
%     zlim([x_min x_max])
%     xlim([y_min y_max])
%     ylim([z_min z_max])
%     view(-55,15);
%     xlabel('Y')
%     ylabel('Z')
%     zlabel('X')
%     subplot (2,3,2)
%     hist(Ecm_step)
%     xlabel('E cm step')
%     subplot (2,3,3)
%     hist(gamma_cm_step)
%     xlabel('\gamma cm step')
%      subplot (2,3,[5 6])
%     hist(Ecm_tot)
%     xlabel('E cm total')
%     filename = [save_dir 'anim_EXP_fig_' num2str(ifig)];
%     fname = [ filename '.png'];
%     print('-dpng', fname);

ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end


set(gca,'FontSize',16)
hold on
%    scatter3(y_min:delta_y:y_max,z_min:delta_z:z_max,x_min:delta_x:x_max,'k')
scatter3(beam_1l(:,6),beam_1l(:,4),beam_1l(:,5),'.b')
scatter3(beam_2l(:,6),beam_2l(:,4),beam_2l(:,5),'.r')
hold off
grid on
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
% view(-55,15);
xlabel('Y')
ylabel('Z')
zlabel('X')
filename = [save_dir '2anim_EXP_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);




% %%
% figure(500)
% hold on
% %    scatter3(y_min:delta_y:y_max,z_min:delta_z:z_max,x_min:delta_x:x_max,'k')
% scatter3(beam_1(:,6),beam_1(:,4),beam_1(:,5),'.b')
% scatter3(beam_2(:,6),beam_2(:,4),beam_2(:,5),'.r')
% hold off
% grid on
% zlim([x_min x_max])
% xlim([y_min y_max])
% ylim([z_min z_max])
% view(-55,15);
% xlabel('Y')
% ylabel('Z')
% zlabel('X')