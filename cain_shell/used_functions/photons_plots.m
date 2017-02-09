function photons_plots(out_folder)
global rflags

if rflags.PLOTS==1
close all;

%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss

photons_data=load([out_folder 'photon_data_plots/photons_data.dat']);%dlmread([out_folder 'cain_tmp/cain_output_photons_' num2str(ni) '.dat'],'',1,0);%read electrons from cain

x_phot=photons_data(:,5);
y_phot=photons_data(:,6);
z_phot=photons_data(:,7);
xp_phot=photons_data(:,9);
yp_phot=photons_data(:,10);
weigth=photons_data(1,3);

full_spectrum=photons_data(:,8)./1e3;
% phot_angle=[phot_angle;sign(photons_data(:,9)).*atan(sqrt(photons_data(:,9).^2+photons_data(:,10).^2)./photons_data(:,11))];
phot_angle=sign(photons_data(:,9)).*atan(sqrt(photons_data(:,9).^2+photons_data(:,10).^2)./photons_data(:,11));

number_of_photons=length(full_spectrum)*weigth;
%SPEED_OF_LIGHT=3e8;

ifig=1;

figure(ifig)
ifig=ifig+1;
hist(full_spectrum,50)
title(['Number photons ' num2str(number_of_photons,'%10.2e')])
filename = [ out_folder 'full_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);

bandwith_cm=[];
num_phot_th=[];
theta_angle=5e-6;
diapason=(10:1:1001);
el_angel=0;
qq=0; num_phot_th=[]; bandwith_cm=[]; bandwith_non_norm=[];
for ni=diapason;
 qq=qq+1;
 el_angel_1=ni*theta_angle;
 aa_1=find(abs(phot_angle)<el_angel_1);
 
 num_phot_th(qq)=length(full_spectrum(aa_1))*weigth;
 bandwith_cm(qq)=std(full_spectrum(aa_1))/mean(full_spectrum(aa_1));
 
 if bandwith_cm(qq)>=0.0048 && bandwith_cm(qq)<=51 % find angle for given bandwith
 el_angel=el_angel_1;
% stop;
 end
end
aa=find(abs(phot_angle)<el_angel);

bandwith=std(full_spectrum(aa))/mean(full_spectrum(aa));
el_angel;
phot_in_bandwidth=length(find(full_spectrum(aa)))*weigth;

%%%%
thetam=el_angel%19.25e-5;
aam=find(abs(phot_angle)<thetam);
bandwith_in_thetam=std(full_spectrum(aam))/mean(full_spectrum(aam));
num_in_thetam=length(find(full_spectrum(aam)))*weigth;


nbin_plot=50;
no_plot_here=1;

xs=linspace(min(full_spectrum(aam)),max(full_spectrum(aam)),nbin_plot);
ys=smooth(hist(full_spectrum(aam),nbin_plot));

xs(find(ys==max(ys)))
pick_of_spectr=xs(find(ys==max(ys)));

if length(pick_of_spectr)>1
pick_of_spectr=mean(pick_of_spectr)
end

std_of_spectr=std(full_spectrum(aam));
mean_of_spectr=mean(full_spectrum(aam));

if no_plot_here==1

figure(ifig)
ifig=ifig+1;
hold on
hist(full_spectrum(aam),nbin_plot,'FaceColor','g')
plot(xs,ys,'-r','LineWidth',3)
plot(xs(find(ys==max(ys))),max(ys),'or','LineWidth',3)
plot(mean_of_spectr,max(ys),'xm','LineWidth',3)
hold off
grid on
set(gca,'FontSize',16)
ylabel('number of scattered macro photons')
xlabel('photons energy (KeV)')
title({['\Theta=' num2str(thetam) ];['  pick of spectr=' num2str(pick_of_spectr,'%1.12e') ];['mean of spectr=' num2str(mean_of_spectr,'%1.12e') ]});
filename = [ out_folder 'pick_photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);
end

%plot(spectr)








std_of_spectr=std(full_spectrum(aam));
mean_of_spectr=mean(full_spectrum(aam));


save([out_folder 'photon_data_plots/phot_spectrum.dat'],'full_spectrum','phot_angle','thetam');

save([out_folder 'photon_data_plots/phot_mat.dat'],'number_of_photons','el_angel','bandwith',...
'phot_in_bandwidth','thetam','bandwith_in_thetam','num_in_thetam','pick_of_spectr','std_of_spectr','mean_of_spectr');


no_plot_here=1;

if no_plot_here==1

xedges = linspace(-1.5e-4,1.5e-4,1e2); yedges =linspace(-1.5e-4,1.5e-4,1e2);
histmat = hist2(x_phot,y_phot, xedges, yedges);

figure(ifig)
ifig=ifig+1;
set(gca,'FontSize',16)
hold on
% mesh(xedges,yedges,histmat')
set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
xlabel('x')
ylabel('y')

filename = [ out_folder 'photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);




figure(ifig)
ifig=ifig+1;
set(gca,'FontSize',16)
subplot(2,1,1)
plot(diapason.*theta_angle,num_phot_th,'-xb')
grid on
set(gca,'FontSize',16)
ylabel('number photons ')
subplot(2,1,2)
hold on
plot(diapason.*theta_angle,bandwith_cm,'--xb')
plot(el_angel,bandwith,'or','LineWidth',2)
hold off
set(gca,'FontSize',16)
ylabel('bandwidth')
xlabel('Theta')
grid on
set(gca,'FontSize',16)
suptitle({['bandwith=' num2str(bandwith) ' for \theta=' num2str(el_angel) ];...
 ['full flux=' num2str(number_of_photons,'%10.2e') ';'...
 ' Flux in \theta =' num2str(length(find(full_spectrum(aa)))*weigth,'%10.2e') ]})
filename = [ out_folder 'photons_plot_' num2str(ifig,'%10.2e') ];
fname = [ filename '.png'];
print('-dpng', fname);





figure(ifig)
ifig=ifig+1;
plot(linspace(0,max(full_spectrum)),histc(full_spectrum(aa),linspace(0,max(full_spectrum)))*weigth,'-r','LineWidth',2)
grid on
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')
% suptitle(['For theta<' num2str(el_angel) ' [rad]']);
suptitle({[''];['For theta<' num2str(el_angel) ' [rad]'];....
    ['Total number of photons=' num2str(number_of_photons,'%10.2e') ];...
    ['number scattered photons in theta<' num2str(length(full_spectrum(aa))*weigth,'%10.2e') ];...
    ['bandwidth=' num2str(std(full_spectrum(aa))/mean(full_spectrum(aa))) ]})
% text(93.1658986175114, 110886076.949367,{['bandwidth=' num2str(std(full_spectrum(aa))/mean(full_spectrum(aa)))];...
%  ['wht']})
% %     ['FWHM=' num2str(0)]})
filename = [ out_folder 'photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);


figure(ifig)
ifig=ifig+1;
subplot(1,2,1)
hold on
hist(phot_angle,50)
hold off
grid on
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('theta [rad]')
% xlabel('photons energy (KeV)')
subplot(1,2,2)
hold on
hist(full_spectrum,50)
hold off
grid on
set(gca,'FontSize',16)
ylabel('number of scattered photons')
%xlabel('theta')
xlabel('photons energy (KeV)')
suptitle({['For all theta'];['Total number of photons=' num2str(number_of_photons,'%10.2e') ]})
filename = [ out_folder 'photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);


figure(ifig)
ifig=ifig+1;
subplot(1,2,1)
hold on
hist(phot_angle(aa),50)
hold off
grid on
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('theta')
% xlabel('photons energy (KeV)')
subplot(1,2,2)
hold on
hist(full_spectrum(aa),50)
hold off
grid on
set(gca,'FontSize',16)
ylabel('number of scattered photons')
%xlabel('theta')
xlabel('photons energy (KeV)')
suptitle({['For theta<' num2str(el_angel) ' [rad]'];['Total number of photons=' num2str(number_of_photons,'%10.2e') ]...
    ;['number scattered photons in theta<' num2str(length(full_spectrum(aa))*weigth,'%10.2e') ]})

filename = [ out_folder 'photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);



nbin_plot=30;
figure(ifig)
ifig=ifig+1;
hold on
plot(linspace(0,max(full_spectrum),nbin_plot),smooth(hist(full_spectrum,nbin_plot)*weigth),'-r','LineWidth',2)
hold off
grid on
ylim([0 max(smooth(hist(full_spectrum,nbin_plot)*weigth))])
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')
title({['For all theta'];['Total number of photons=' num2str(number_of_photons,'%10.2e') ]})
filename = [ out_folder 'photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);
end
%%%%%%
% 
% figure(ifig)
% ifig=ifig+1;
% hold on
% plot(phot_angle,full_spectrum,'.b')
% hold off
% grid on
% set(gca,'FontSize',16)
% ylabel('energy of scattered photons')
% xlabel('scattered angle')
% suptitle({['For all theta'];['Total number of photons=' num2str(number_of_photons,'%10.2e') ]});
% filename = [ out_folder 'photons_plot_' num2str(ifig) ];
% fname = [ filename '.png'];
% print('-dpng', fname);



% 
% figure(ifig)
% ifig=ifig+1;
% hold on
% plot(phot_angle(aa),full_spectrum(aa),'.b')
% hold off
% grid on
% set(gca,'FontSize',16)
% ylabel('energy of scattered photons')
% xlabel('scattered angle')
% suptitle({['For theta<' num2str(el_angel) ' [rad]'];['Total number of photons=' num2str(number_of_photons,'%10.2e') ];['number scattered photons in theta<' num2str(length(full_spectrum(aa))*weigth,'%10.2e') ]})
% filename = [ out_folder 'photons_plot_' num2str(ifig) ];
% fname = [ filename '.png'];
% print('-dpng', fname);

% 
% 
% 
% ifig=ifig+1;
% figure(ifig)
% subplot 221
% set(gca,'FontSize',16)
% plot(x_phot,y_phot,'.b')
% xlabel('x [m]')
% ylabel('y [m]')
% subplot 222
% set(gca,'FontSize',16)
% plot(x_phot,xp_phot,'.b')
% xlabel('x [m]')
% ylabel('Px [eV/c]')
% subplot 223
% set(gca,'FontSize',16)
% plot(y_phot,yp_phot,'.b')
% xlabel('y [m]')
% ylabel('Px [eV/c]')
% subplot 224
% set(gca,'FontSize',16)
% plot(xp_phot,yp_phot,'.b')
% xlabel('Px [eV/c]')
% ylabel('Py [eV/c]')
% suptitle({['For all theta'];['Total number of photons=' num2str(number_of_photons,'%10.2e') ]});
% filename = [ out_folder 'photons_plot_' num2str(ifig) ];
% fname = [ filename '.png'];
% print('-dpng', fname);


% 
% ifig=ifig+1;
% figure(ifig)
% 
% subplot 221
% set(gca,'FontSize',16)
% plot(x_phot(aa),y_phot(aa),'.b')
% xlabel('x [m]')
% ylabel('y [m]')
% subplot 222
% set(gca,'FontSize',16)
% plot(x_phot(aa),xp_phot(aa),'.b')
% xlabel('x [m]')
% ylabel('Px [eV/c]')
% subplot 223
% set(gca,'FontSize',16)
% plot(y_phot(aa),yp_phot(aa),'.b')
% xlabel('y [m]')
% ylabel('Px [eV/c]')
% subplot 224
% set(gca,'FontSize',16)
% plot(xp_phot(aa),yp_phot(aa),'.b')
% xlabel('Px [eV/c]')
% ylabel('Py [eV/c]')
% suptitle({['For theta<' num2str(el_angel) ' [rad]'];['Total number of photons=' num2str(number_of_photons,'%10.2e') ];['number scattered photons in theta<' num2str(length(full_spectrum(aa))*weigth,'%10.2e') ]})
% filename = [ out_folder 'photons_plot_' num2str(ifig) ];
% fname = [ filename '.png'];
% print('-dpng', fname);
% 
% 
% 
% ifig=ifig+1;
% figure(ifig)
% subplot 221
% set(gca,'FontSize',16)
% hold on
% plot(x_phot,y_phot,'.b')
% plot(x_phot(aa),y_phot(aa),'.r')
% hold off
% xlabel('x [m]')
% ylabel('y [m]')
% subplot 222
% hold on
% set(gca,'FontSize',16)
% plot(x_phot,xp_phot,'.b')
% plot(x_phot(aa),xp_phot(aa),'.r')
% hold off
% xlabel('x [m]')
% ylabel('Px [eV/c]')
% subplot 223
% hold on
% set(gca,'FontSize',16)
% plot(y_phot,yp_phot,'.b')
% plot(y_phot(aa),yp_phot(aa),'.r')
% hold off
% xlabel('y [m]')
% ylabel('Px [eV/c]')
% subplot 224
% set(gca,'FontSize',16)
% hold on
% plot(xp_phot,yp_phot,'.b')
% plot(xp_phot(aa),yp_phot(aa),'.r')
% 
% hold off
% xlabel('Px [eV/c]')
% ylabel('Py [eV/c]')
% 
% suptitle({['For theta<' num2str(el_angel) ' [rad]'];['Total number of photons=' num2str(number_of_photons,'%10.2e') ];['number scattered photons in theta<' num2str(length(full_spectrum(aa))*weigth,'%10.2e') ]})
% 
% filename = [ out_folder 'photons_plot_' num2str(ifig) ];
% fname = [ filename '.png'];
% print('-dpng', fname);



%%%%%%%%%%%%%



%nbin=[101];
%% xedges = linspace(-2,2,nbin); yedges = linspace(0,50,nbin);
%xedges = linspace(-1,1,nbin); yedges = linspace(0,max(full_spectrum/1e3),nbin);
%histmat = hist2(phot_angle*1e3, full_spectrum/1e3, xedges, yedges);
%
%new_m=zeros(nbin,nbin);
%for nni=1:1:nbin-1
%    new_m(nni,:)=histmat(nni,:)./((pi/2).*(abs((xedges(nni+1)).^2-(xedges(nni)).^2)));
%    
%end
%ifig=ifig+1;
%figure(ifig)
%
%set(pcolor(xedges,yedges,new_m'),'EdgeColor','none')
%%  surf(xedges,yedges,new_m')
%% colormap([1 1 1;0.857142865657806 0.857142865657806 1;0.714285731315613 0.714285731315613 1;0.571428596973419 0.571428596973419 1;0.428571432828903 %0.428571432828903 1;0.28571429848671 0.28571429848671 1;0.142857149243355 0.142857149243355 1;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 0.9375;0.125 1 0.875;0.1875 1 0.8125;0.25 1 0.75;0.3125 1 0.6875;0.375 1 0.625;0.4375 1 0.5625;0.5 1 0.5;0.5625 1 0.4375;0.625 1 0.375;0.6875 1 0.3125;0.75 1 0.25;0.8125 1 0.1875;0.875 1 0.125;0.9375 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0;0.5 0 0]);
%shading interp
%colorbar
%set(gca,'FontSize',16)
%title({['']; })
%set(gca,'FontSize',16)
%xlabel('angle (mrad)');
%ylabel('photons energy (MeV)');
%filename = [ out_folder 'photons_plot_' num2str(ifig) ];
%fname = [ filename '.png'];
%print('-dpng', fname);




end