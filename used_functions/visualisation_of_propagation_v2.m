function visualisation_of_propagation_v2(beam_1l,beam_2l,x_max,x_min,y_max,y_min,z_max,z_min,delta_x,delta_y,delta_z)
global ifig save_dir rflags
if rflags.PLOTS ==1;
    figure(ifig)
    ifig=ifig+1;
set(gca,'FontSize',16)
hold on
%    scatter3(y_min:delta_y:y_max,z_min:delta_z:z_max,x_min:delta_x:x_max,'k')
scatter3(beam_1l(:,6),beam_1l(:,4),beam_1l(:,5),'.b')
scatter3(beam_2l(:,6),beam_2l(:,4),beam_2l(:,5),'.r')
% scatter3(y_max,z_max,x_max,'.r')
% scatter3(y_min,z_min,x_min,'.b')
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
 view(-40,90);
xlabel('Y')
ylabel('Z')
zlabel('X')
    
    % %     subplot(121)
%     hold on
%     %    scatter3(y_min:delta_y:y_max,z_min:delta_z:z_max,x_min:delta_x:x_max,'k')
%     scatter3(beam_1l(:,6),beam_1l(:,4),beam_1l(:,5),'.b')
%     scatter3(beam_2l(:,6),beam_2l(:,4),beam_2l(:,5),'.r')
%     hold off
%     grid on
% %     ax = gca;
% %     ax.XTick = [y_min:delta_y:y_max];
% %     ay = gca;
% %     ay.YTick = [z_min:delta_z:z_max];
% %     az = gca;
% %     az.ZTick = [x_min:delta_x:x_max];
%     zlim([x_min x_max])
%     xlim([y_min y_max])
%     ylim([z_min z_max])
% %     view(-90,90);
%     xlabel('Y')
%     ylabel('Z')
%     zlabel('X')
% %     subplot(122)
% %     hist(Ecm_step)
% %     xlabel('E cm')
    filename = [save_dir 'A_EXP_fig_' num2str(ifig)];
    fname = [ filename '.png'];
    print('-dpng', fname);
end