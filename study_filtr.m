clear all; close all; clc;
rflags.PLOTS =1;
ifig=0;
% stop
save_dir=[pwd '/GGfel_low_10k_3sigma_FEL_GG_cut_1/nbin_x_21_y_21_z_21/'];
filename=[save_dir 'Num_events_'];
load([save_dir 'main_out_put.dat'],'-mat');
ad=load([save_dir 'out_put.dat'],'-mat');
% pair_info=[Vq E_3 E_4 theta_3 theta_4 phi_3 phi_4 cos_alpha cos_theta Ecm_pair gamma_cm_pair];
%%

sc2=load('/Users/drilrrr/Documents/ROSE_CODE_production/FEL_COMPTON_BG_17_01_2016/BG_FEL_compton_table_low_10k_3sigma_FEL_GG_cut_1/nbin_x_21_y_21_z_21/main_out_put.dat','-mat');


ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
histogram2(sc2.pair_full(:,3),sc2.pair_full(:,5),[150 150])
xlabel('Energy','FontSize',20)
ylabel('\theta','FontSize',20)
zlabel('# \gamma','FontSize',20)
%
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);





%%
photons_data=dlmread(['/Users/drilrrr/Documents/ROSE_CODE_production/GAMMA_GAMMA/Channeling_2016/CAIN_GG/BIG_STAT_check_10K_no_halo_sigma_3_defocus_1/photon_data_plots/photons_data.dat'],'',0,0);%read photons from cain
% photons_data=dlmread([pwd '/gamma_gamma_low/photons_data_10K_3sigma.dat'],'',0,0);%read photons from cain
full_spectrum=photons_data(:,8)./1e3;
% phot_angle=[phot_angle;sign(photons_data(:,9)).*atan(sqrt(photons_data(:,9).^2+photons_data(:,10).^2)./photons_data(:,11))];
phot_angle=abs(atan(sqrt(photons_data(:,9).^2+photons_data(:,10).^2)./photons_data(:,11)));
%% histogram beam1 vs theta and Energy_lab
ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
hold on
pr_com=histogram2(phot_angle,full_spectrum,[20 20])
hold off
ylabel('Energy','FontSize',20)
xlabel('\theta','FontSize',20)
zlabel('# \gamma_1 ','FontSize',20)
%
pr_com_bin=pr_com.Values.*photons_data(1,3);
X_c_com=linspace(min(full_spectrum),max(full_spectrum),20);
Y_c_com=linspace(min(phot_angle),max(phot_angle),20);

ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
surf(X_c_com,Y_c_com,pr_com_bin')

%%

%% histogram beam1 vs theta and Energy_lab
ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
histogram2(pair_full(:,2),pair_full(:,4),[50 50])
xlabel('Energy','FontSize',20)
ylabel('\theta','FontSize',20)
zlabel('# \gamma_3 ','FontSize',20)
%
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);



%%
theta_cut=pi-0.05;
fan1=find(pair_full(:,4)<theta_cut);
fan2=find(pair_full(:,4)>theta_cut);
ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
subplot 121
histogram2(pair_full(fan1,2),pair_full(fan1,4),[50 50])
%     histogram2(pair_full(fan1,3),pair_full(fan1,5),[50 50])
xlabel('Energy','FontSize',20)
ylabel('\theta','FontSize',20)
zlabel('# \gamma_3 ','FontSize',20)
title(['for \theta<' num2str(theta_cut)])
subplot 122
histogram2(pair_full(fan2,2),pair_full(fan2,4),[50 50])
%     histogram2(pair_full(fan1,3),pair_full(fan1,5),[50 50])
xlabel('Energy','FontSize',20)
ylabel('\theta','FontSize',20)
zlabel('# \gamma_3 ','FontSize',20)
title(['for \theta>' num2str(theta_cut)])


ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
subplot 121
histogram2(pair_full(fan1,3),pair_full(fan1,5),[50 50])
xlabel('Energy','FontSize',20)
ylabel('\theta','FontSize',20)
zlabel('# \gamma_4 ','FontSize',20)
title(['for \theta<' num2str(theta_cut)])
subplot 122
histogram2(pair_full(fan2,3),pair_full(fan2,5),[50 50])
xlabel('Energy','FontSize',20)
ylabel('\theta','FontSize',20)
zlabel('# \gamma_4 ','FontSize',20)
title(['for \theta>' num2str(theta_cut)])

%  stop
%% histogram beam 2 vs theta and Energy_lab
ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
histogram2(pair_full(:,3),pair_full(:,5),[50 50])
xlabel('Energy','FontSize',20)
ylabel('\theta','FontSize',20)
zlabel('# \gamma_2 ','FontSize',20)
%
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);




%% histogram N pair vs theta and Energy_lab
ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
as=histogram2([pair_full(:,2); pair_full(:,3)],[pair_full(:,4); pair_full(:,5)],[50 50])
xlabel('Energy','FontSize',20)
ylabel('\theta','FontSize',20)
zlabel('# of pair','FontSize',20)
%
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);

N_numer=sum(sum(as.Values,2));

%%
n_bin=20;

%% Make a histogram of E_CoM

ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
hist(ad.Ecm_tot,n_bin)
xlabel('E_{CoM}','FontSize',20)
ylabel('Number of pair','FontSize',20)
%% and extract from it coordinates of center of bin

[N_pair_per_bin cc_bin]=hist(ad.Ecm_tot,n_bin);

N_gammas=2*N_pair_per_bin*ad.weight_1*ad.weight_2;

ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
bar(cc_bin,N_gammas)
xlabel('E_{CoM}','FontSize',20)
ylabel('Number of gammas','FontSize',20)
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);


%% load diff crossection
ev=load([pwd '/used_functions/rs_grid.dat'],'ascii')';
cs=load([pwd '/used_functions/w_grid.dat'],'ascii')';
cross=load([pwd '/used_functions/cs_grid.dat'],'ascii')';



ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
mesh((cs(:,1)),ev(1,:),cross')
view(90, 0)
zlabel('cross-section','FontSize',20)
ylabel('E_{CoM}','FontSize',20)
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);


% ifig=ifig+1;
% figure(ifig);
% plot(ev(1,:),(2*pi/4*pi)*sum(cross(:,1:1001)).*(cs(101,1)-cs(100,1)),'LineWidth',2)
% ylabel('cross-section','FontSize',20)
% xlabel('E_{CoM}','FontSize',20)
% fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
% print('-r300','-dpng', fname2); %print('-depsc', fname);


tot_e_cross=(pi)*sum(cross(:,1:1001)).*(cs(101,1)-cs(100,1));
en_x=ev(1,:);

ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
plot(ev(1,:),tot_e_cross,'LineWidth',2)
ylabel('cross-section','FontSize',20)
xlabel('E_{CoM}','FontSize',20)
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);


% ev_l=linspace(0.1,2.4,n_bin);
ev_l=cc_bin/1e6;

% ev_bin = interp1(ev_l,1:numel(ev_l),tot_e_cross,'nearest')';
VCE=interpn(en_x,tot_e_cross,ev_l,'cubic');

ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
hold on
bar(ev_l,VCE)
plot(en_x,tot_e_cross,'-',ev_l,VCE,'o')
hold off
ylabel('cross-section','FontSize',20)
xlabel('E_{CoM}','FontSize',20)
title('bar and dots are interpolated data')
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);

%%
ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
% Create figure
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
% Create bar
bar(ev_l,VCE);
% Create loglog
loglog(en_x,tot_e_cross);
% Create loglog
loglog(ev_l,VCE,'Marker','o','LineStyle','none');
% Create xlabel
xlabel('E_{CoM}');
% Create ylabel
ylabel('cross-section');
% Create title
title('bar and dots are interpolated data');
% Set the remaining axes properties
set(axes1,'XScale','log','YScale','log');
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);
%%
ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
bar(ev_l,(VCE.*1e-34./(ad.delta_x.*ad.delta_y)).*N_gammas)

ylabel('N event','FontSize',20)
xlabel('E_{CoM}','FontSize',20)
title(['sum ' num2str(sum((VCE.*1e-34./(ad.delta_x.*ad.delta_y)).*N_gammas)) ])
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);




N_tot=sum((VCE.*1e-34./(ad.delta_x.*ad.delta_y)).*N_gammas);

Aaa=N_tot/N_numer;

real_bin=as.Values.*Aaa;
X_c=linspace(min([pair_full(:,2); pair_full(:,3)]),max([pair_full(:,2); pair_full(:,3)]),50);
Y_c=linspace(min([pair_full(:,4); pair_full(:,5)]),max([pair_full(:,4); pair_full(:,5)]),50);


%% histogram N pair vs theta and Energy_lab
ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
surf(X_c,Y_c,real_bin')
%=histogram2([pair_full(:,2); pair_full(:,3)],[pair_full(:,4); pair_full(:,5)],[5 5])
xlabel('Energy','FontSize',20)
ylabel('\theta','FontSize',20)
zlabel('Nubmer of \gamma','FontSize',20)
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);




ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
hold on
%  colormap([1 0 0;0 0 1]) %red and blue
surf(X_c,Y_c,(as.Values./max(max(as.Values)))','EdgeColor','none')
surf(X_c_com.*1e3,pi-Y_c_com,(pr_com.Values./max(max(pr_com.Values))),'EdgeColor','none')
hold off


% sc2en_all=sc2.pair_full(:,3);
% sc2theta_all=sc2.pair_full(:,5);

small_energy=find(sc2.pair_full(:,3)<1.2e6);

sc2en=sc2.pair_full(small_energy,3);
sc2theta=sc2.pair_full(small_energy,5);

ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
sc2h=histogram2(sc2en,pi-sc2theta,[50 50])
xlabel('Energy','FontSize',20)
ylabel('\theta','FontSize',20)
zlabel('# \gamma','FontSize',20)
%
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); %print('-depsc', fname);

X_c_sc2=linspace(min(sc2en),max(sc2en),50);
Y_c_sc2=linspace(min(sc2theta),max(sc2theta),50);


ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
hold on
%  colormap([1 0 0;0 0 1]) %red and blue
 surf(X_c,Y_c,(as.Values./max(max(as.Values)))','EdgeColor','none')
 surf(X_c_com.*1e3,pi-Y_c_com,(pr_com.Values./max(max(pr_com.Values))),'EdgeColor','none')
 surf(X_c_sc2,Y_c_sc2,(sc2h.Values./max(max(sc2h.Values)))','EdgeColor','none')
hold off
xlabel('Energy','FontSize',20)
ylabel('\theta','FontSize',20)



