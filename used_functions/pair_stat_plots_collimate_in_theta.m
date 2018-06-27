function pair_stat_plots_collimate_in_theta(pair_full,theta_max)

global ifig save_dir
global rflags 
filename=[save_dir 'pair_info_'];
if rflags.e_gamma==0


col=find(pair_full(:,5)<theta_max);
theta=pair_full(col,5);
en=pair_full(col,3);


ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
as=histogram2(pair_full(col,3),pair_full(col,5),[50 50],'DisplayStyle','tile','ShowEmptyBins','on');
xlabel('Energy','FontSize',20)
ylabel('\theta','FontSize',20)
zlabel('# \gamma','FontSize',20)
view(90,-90)
title(['in \theta_{max}=' num2str(theta_max)])
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);
%

ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
histogram(pair_full(col,3));
xlabel('Energy','FontSize',20)
title(['in \theta_{max}=' num2str(theta_max)])
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); 


%% 
nbin=50;
xedges = linspace(0,max(theta),nbin); yedges =linspace(min(en),max(en),nbin);
histmat = hist2(theta,en, xedges, yedges);

ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
set(gca,'FontSize',16)
% mesh(xedges,yedges,histmat')
set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
ylabel('Energy','FontSize',20)
xlabel('\theta','FontSize',20)
title(['in \theta_{max}=' num2str(theta_max)])


end

