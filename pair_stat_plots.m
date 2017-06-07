function pair_stat_plots(pair_full)

global ifig save_dir
global rflags 
filename=[save_dir 'pair_info_'];
if rflags.e_gamma==0


if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
ifig=ifig+1;
set(gca,'FontSize',12)
hist(pair_full(:,3),20)
xlabel('E_{CoM}','FontSize',20)
ylabel('Scattered photons','FontSize',20)
% title(['N \gamma  ' num2str(3.52e07,'%10.2e') ])
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);


if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
ifig=ifig+1;
set(gca,'FontSize',12)
hist(pair_full(:,5),20)
xlabel('\theta','FontSize',20)
ylabel('Scattered photons','FontSize',20)
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);


%% histogram beam 2 vs theta and Energy_lab
ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
as=histogram2(pair_full(:,3),pair_full(:,5),[50 50],'DisplayStyle','tile','ShowEmptyBins','on');
xlabel('Energy','FontSize',20)
ylabel('\theta','FontSize',20)
zlabel('# \gamma','FontSize',20)
view(90,-90)
%
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);

%% histogram beam 2 vs theta and Energy_lab

indc=find(pair_full(:,5)<5);

ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
as=histogram2(pair_full(indc,3),pair_full(indc,5),[50 50],'DisplayStyle','tile','ShowEmptyBins','on');
xlabel('Energy','FontSize',20)
ylabel('\theta','FontSize',20)
zlabel('# \gamma','FontSize',20)
view(90,-90)
%
% title(['N \gamma  ' num2str(3.52e07,'%10.2e') ])
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);



end

