function Plot_for_cheking_sizes_after_l_max(beam_1,beam_2,l_max)

%% Plot for cheking sizes after l_max
global rflags ifig save_dir_start

if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
ifig=ifig+1;
subplot 221
plot(beam_1(:,5),beam_1(:,6),'.r')
grid on
subplot 222
plot(beam_2(:,5),beam_2(:,6),'.r')
grid on
[beam_1l]=beam_drift(beam_1,l_max);
[beam_2l]=beam_drift(beam_2,-l_max);
subplot 223
plot(beam_1l(:,5),beam_1l(:,6),'.b')
grid on
subplot 224
plot(beam_2l(:,5),beam_2l(:,6),'.b')
grid on
filename = [save_dir_start 'gamma_EXP_fig_' num2str(ifig)];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
% print('-depsc', fname);
print('-r300','-dpng', fname2);



%%
xedges_1=[]; yedges_1=[];
xedges_1 = linspace(min(beam_1(:,5)),max(beam_1(:,5)),1e2); yedges_1 =linspace(min(beam_1(:,6)),max(beam_1(:,6)),1e2);
histmat_1 = hist2(beam_1(:,5),beam_1(:,6), xedges_1, yedges_1);
xedges_2=[]; yedges_2=[];
xedges_2 = linspace(min(beam_2(:,5)),max(beam_2(:,5)),1e2); yedges_2 =linspace(min(beam_2(:,6)),max(beam_2(:,6)),1e2);
histmat_2 = hist2(beam_2(:,5),beam_2(:,6), xedges_2, yedges_2);
xedges_1l=[]; yedges_1l=[];
xedges_1l = linspace(min(beam_1l(:,5)),max(beam_1l(:,5)),1e2); yedges_1l =linspace(min(beam_1l(:,6)),max(beam_1l(:,6)),1e2);
histmat_1l = hist2(beam_1l(:,5),beam_1l(:,6), xedges_1l, yedges_1l);
xedges_2l=[]; yedges_2l=[];
xedges_2l = linspace(min(beam_2l(:,5)),max(beam_2l(:,5)),1e2); yedges_2l =linspace(min(beam_2l(:,6)),max(beam_2l(:,6)),1e2);
histmat_2l = hist2(beam_2l(:,5),beam_2l(:,6), xedges_2l, yedges_2l);







if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
ifig=ifig+1;
subplot 221
imagesc(xedges_1,yedges_1,histmat_1')%,'EdgeColor','none')
grid on
title('beam 1')
subplot 222
imagesc(xedges_2,yedges_2,histmat_2')%,'EdgeColor','none')
grid on
title('beam 2')
% [beam_1l]=beam_drift(beam_1,l_max);
% [beam_2l]=beam_drift(beam_2,-l_max);
subplot 223
imagesc(xedges_1l,yedges_1l,histmat_1l')%,'EdgeColor','none')
grid on
title('beam 1+l')
subplot 224
imagesc(xedges_2l,yedges_2l,histmat_2l')%,'EdgeColor','none')
grid on
title('beam 2+l')
filename = [save_dir_start 'gamma_EXP_fig_' num2str(ifig)];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
% print('-depsc', fname);
print('-r300','-dpng', fname2);

% stop


