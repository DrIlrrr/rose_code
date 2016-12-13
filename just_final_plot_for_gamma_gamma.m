clear all; close all; clc;
rflags.PLOTS =1;
ifig=0;
stop
save_dir=[pwd '/GG_test_0.001_GG2/nbin_x_21_y_21_z_21/'];
filename=[save_dir 'Num_events_'];
load([save_dir 'main_out_put.dat'],'-mat');
ad=load([save_dir 'out_put.dat'],'-mat');
% pair_info=[Vq E_3 E_4 theta_3 theta_4 phi_3 phi_4 cos_alpha cos_theta Ecm_pair gamma_cm_pair];

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

%% histogram N pair Energy_lab
ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
histogram([pair_full(:,2); pair_full(:,3)],20)
xlabel('Energy','FontSize',20)
% ylabel('\theta','FontSize',20)
ylabel('# of pair','FontSize',20)
%
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); print('-depsc', fname);

%% histogram N pair vs theta
ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
histogram([pair_full(:,4); pair_full(:,5)],50)
% xlabel('Energy','FontSize',20)
xlabel('\theta','FontSize',20)
ylabel('# of pair','FontSize',20)
%
fname = [ filename num2str(ifig) '.eps'];fname2 = [ filename num2str(ifig) '.png'];
print('-r300','-dpng', fname2); print('-depsc', fname);

stop
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
ev=load('/Users/drilrrr/Dropbox/rabota/rose_code/used_functions/rs_grid.dat','ascii')';
cs=load('/Users/drilrrr/Dropbox/rabota/rose_code/used_functions/w_grid.dat','ascii')';
cross=load('/Users/drilrrr/Dropbox/rabota/rose_code/used_functions/cs_grid.dat','ascii')';



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

