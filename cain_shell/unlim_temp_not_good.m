% one beam several events
clear all; close all; clc;
stop
make_path
start_date=datestr(now);
%for new global
global home_dir DIRECTORY_FOR_CAIN BASE_DIRECTORY;
global rflags
[rflags] = flags_for_run;
home_dir=[pwd '/CAIN/'];

ifig=1;
rflags.PLOTS =1;
just_plots=1;
qq=0;
scan_spectrum=[];
phot_angle=[];


STOKES(1,:)=[1 0 0];
STOKES(2,:)=[0 1 0];
STOKES(3,:)=[0 0 1];

aaa=[1 2 3]
for var_for_scan=aaa;
    qq=qq+1;
    
    
    rflags.STOKES=STOKES(var_for_scan,:);
    
    
    
    
    rflags.pulseE=1; %laser puse energy [J]
    
    
    
    defoc_param=1.5;%defocusing parameter by defoult 1 is no defocusing
    BASE_DIRECTORY = [pwd '/dif_STOKES_' num2str(rflags.STOKES(1)) '_' num2str(rflags.STOKES(2)) '_' num2str(rflags.STOKES(3)) '_polar_pulseE_' num2str(rflags.pulseE) '_hi_600_defocus_' num2str(defoc_param) '_PulseE_400mJ/'];
    
    [full_spectrum,phot_angle,weigth,x_phot,y_phot,z_phot,xp_phot,yp_phot,zp_phot,Sx,Sy,Sz]=scan_photons_plots(BASE_DIRECTORY);
    
    eval(['full_spectrum_' int2str(qq) '=full_spectrum;'])
    eval(['phot_angle_' int2str(qq) '=phot_angle;'])
    eval(['x_coor' int2str(qq) '=x_phot;'])
    eval(['y_coor' int2str(qq) '=y_phot;'])
    eval(['z_coor' int2str(qq) '=z_phot;'])
    eval(['xp_' int2str(qq) '=xp_phot;'])
    eval(['yp_' int2str(qq) '=yp_phot;'])
    eval(['zp_' int2str(qq) '=zp_phot;'])
    eval(['Sx_' int2str(qq) '=Sx;'])
    eval(['Sy_' int2str(qq) '=Sy;'])
    eval(['Sz_' int2str(qq) '=Sz;'])
end




ww=0;
for ni=aaa
    ww=ww+1
    eval(['aa_' num2str(ww) '=find(abs(phot_angle_' num2str(ww) '));']);
    num_in_b(ww)=eval(['length(full_spectrum_' int2str(ww) '(aa_' int2str(ww) '));']);
    
end

figure(ifig)
ifig=ifig+1;
hold on
plot(aaa,num_in_b*weigth,'-.xb','LineWidth',0.5)
plot(aaa,num_in_b*weigth,'xb','LineWidth',3)
% plot(aaa,NUM_ph,'or','LineWidth',3)
hold off
grid on
ylim([0 max(num_in_b*weigth)+max(num_in_b*weigth)*1e-1])
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('\alpha_0')
% legend('\delta t=0 [ps]','\delta t=1 [ps]','\delta t=2 [ps]',0)
filename = ['EXP_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);






theta=5e-5;
ww=0;
for ni=aaa
    ww=ww+1
    eval(['aa_' num2str(ww) '=find(abs(phot_angle_' num2str(ww) ')<' num2str(theta) ');']);
    num_in_b(ww)=eval(['length(full_spectrum_' int2str(ww) '(aa_' int2str(ww) '));']);
    
end

nbin=50;
min_val_linspase=12000
figure(ifig)
ifig=ifig+1;
hold on
cc=jet(ww)
for ni=1:1:ww
    plot(linspace(min_val_linspase,max(eval(['full_spectrum_' int2str(ni) ])),nbin),...
        histc(eval(['full_spectrum_' int2str(ni) '(aa_' int2str(ni) ')']),...
        linspace(min_val_linspase,max(eval(['full_spectrum_' int2str(ni) ])),nbin))*weigth,'color',cc(ni,:),'LineWidth',2)
end
grid on
set(gca,'FontSize',16)
title({['Spectrum of scattered photons in Theta=' num2str(theta) ' [rad]'],...
    ['bandwith =' num2str(std(full_spectrum_1(aa_1))/mean(full_spectrum_1(aa_1))) '']})
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')
% legend(['\alpha_0=0 ; photons = ' num2str(num_in_b(1)*weigth,'%2.1e') ],...
%     ['\alpha_0=2.5; photons = ' num2str(num_in_b(2)*weigth,'%2.1e') ],...
%     ['\alpha_0=5; photons = ' num2str(num_in_b(3)*weigth,'%2.1e') ],...
%     ['\alpha_0=7.5; photons = ' num2str(num_in_b(4)*weigth,'%2.1e') ],0)
filename = ['EXP_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);




figure(ifig)
ifig=ifig+1;

for ni=1:1:3
    
    xedges=[];
    num_bin=50;
    xedges = linspace(min(eval(['full_spectrum_' int2str(ni) ])),max(eval(['full_spectrum_' int2str(ni) ])),num_bin);...
        yedges =linspace(min(eval(['Sx_' int2str(ni) ])),1,num_bin);
    histmat = hist2(eval(['full_spectrum_' int2str(ni) ]),eval(['Sx_' int2str(ni) ]), xedges, yedges);
    
    subplot(3,3,ni)
    set(gca,'FontSize',16)
    hold on
%     mesh(xedges,yedges,histmat')
    set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    set(gca,'FontSize',16)
    xlabel('Photon energy')
    ylabel('Sx')
    
    
    xedges=[];
    num_bin=50;
    xedges = linspace(min(eval(['full_spectrum_' int2str(ni) ])),max(eval(['full_spectrum_' int2str(ni) ])),num_bin);...
        yedges =linspace(min(eval(['Sy_' int2str(ni) ])),1,num_bin);
    histmat = hist2(eval(['full_spectrum_' int2str(ni) ]),eval(['Sy_' int2str(ni) ]), xedges, yedges);
    
    subplot(3,3,ni+3)
    set(gca,'FontSize',16)
    hold on
%     mesh(xedges,yedges,histmat')
    set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    set(gca,'FontSize',16)
    xlabel('Photon energy')
    ylabel('Sy')
    
    
    xedges=[];
    num_bin=50;
    xedges = linspace(min(eval(['full_spectrum_' int2str(ni) ])),max(eval(['full_spectrum_' int2str(ni) ])),num_bin);...
        yedges =linspace(min(eval(['Sz_' int2str(ni) ])),1,num_bin);
    histmat = hist2(eval(['full_spectrum_' int2str(ni) ]),eval(['Sz_' int2str(ni) ]), xedges, yedges);
    
    subplot(3,3,ni+6)
    set(gca,'FontSize',16)
    hold on
%     mesh(xedges,yedges,histmat')
    set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    set(gca,'FontSize',16)
    xlabel('Photon energy')
    ylabel('Sz')
    
end
filename = ['EXP_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);




figure(ifig)
ifig=ifig+1;

for ni=1:1:3
    
    xedges=[];
    num_bin=20;
    xedges = linspace(min(eval(['full_spectrum_' int2str(ni) ])),max(eval(['full_spectrum_' int2str(ni) ])),num_bin);...
        yedges =linspace(min(eval(['Sx_' int2str(ni) ])),1,num_bin);
    % xedges = linspace(min(eval(['full_spectrum_' int2str(ni)  '(aa_' int2str(ni) ')'])),max(eval(['full_spectrum_' int2str(ni)  '(aa_' int2str(ni) ')'])),num_bin);...
    %     yedges =linspace(min(eval(['Sx_' int2str(ni)  '(aa_' int2str(ni) ')'])),1,num_bin);
    histmat = hist2(eval(['full_spectrum_' int2str(ni)  '(aa_' int2str(ni) ')']),eval(['Sx_' int2str(ni)  '(aa_' int2str(ni) ')']), xedges, yedges);
    
    subplot(3,3,ni)
    set(gca,'FontSize',16)
    hold on
%     mesh(xedges,yedges,histmat')
    set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    set(gca,'FontSize',16)
    % xlabel('Photon energy')
    ylabel('Sx')
    
    
    xedges=[];
    % num_bin=50;
    xedges = linspace(min(eval(['full_spectrum_' int2str(ni) ])),max(eval(['full_spectrum_' int2str(ni) ])),num_bin);...
        yedges =linspace(min(eval(['Sy_' int2str(ni) ])),1,num_bin);
    % xedges = linspace(min(eval(['full_spectrum_' int2str(ni)  '(aa_' int2str(ni) ')'])),max(eval(['full_spectrum_' int2str(ni)  '(aa_' int2str(ni) ')'])),num_bin);...
    %     yedges =linspace(min(eval(['Sy_' int2str(ni)  '(aa_' int2str(ni) ')'])),1,num_bin);
    histmat = hist2(eval(['full_spectrum_' int2str(ni)  '(aa_' int2str(ni) ')']),eval(['Sy_' int2str(ni)  '(aa_' int2str(ni) ')']), xedges, yedges);
    
    subplot(3,3,ni+3)
    set(gca,'FontSize',16)
    hold on
%     mesh(xedges,yedges,histmat')
    set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    set(gca,'FontSize',16)
    % xlabel('Photon energy')
    ylabel('Sy')
    
    
    xedges=[];
    % num_bin=50;
    xedges = linspace(min(eval(['full_spectrum_' int2str(ni) ])),max(eval(['full_spectrum_' int2str(ni) ])),num_bin);...
        yedges =linspace(min(eval(['Sz_' int2str(ni) ])),1,num_bin);
    % xedges = linspace(min(eval(['full_spectrum_' int2str(ni)  '(aa_' int2str(ni) ')'])),max(eval(['full_spectrum_' int2str(ni)  '(aa_' int2str(ni) ')'])),num_bin);...
    %     yedges =linspace(min(eval(['Sz_' int2str(ni)  '(aa_' int2str(ni) ')'])),1,num_bin);
    histmat = hist2(eval(['full_spectrum_' int2str(ni)  '(aa_' int2str(ni) ')']),eval(['Sz_' int2str(ni)  '(aa_' int2str(ni) ')']), xedges, yedges);
    
    subplot(3,3,ni+6)
    set(gca,'FontSize',16)
    hold on
%     mesh(xedges,yedges,histmat')
    set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    set(gca,'FontSize',16)
    xlabel('Photon energy')
    ylabel('Sz')
    
end
suptitle(['in Theta<' num2str(theta) ' [rad]'])
filename = ['EXP_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);



%% momentum plot

figure(ifig)
ifig=ifig+1;
theta=5e-5;
ww=0;
for ni=aaa
    ww=ww+1
    eval(['aaa_' num2str(ww) '=find(abs(phot_angle_' num2str(ww) ')<' num2str(theta) ');']);
    %     num_in_b(ww)=eval(['length(full_spectrum_' int2str(ww) '(aa_' int2str(ww) '));']);
    
end
suptitle(['theta<' num2str(theta)...
    '; bandwith =' num2str(std(full_spectrum_1(aaa_1))/mean(full_spectrum_1(aaa_1)))  ...
    ])
for ni=1:1:3
    xedges=[]; yedges=[]; histmat=[];
    xedges = linspace(-0.6e4,0.6e4,1e2); yedges =linspace(-0.6e4,0.6e4,1e2);
    histmat = hist2(eval(['xp_' int2str(ni) ]),eval(['yp_' int2str(ni) ]), xedges, yedges);
    
    subplot(2,3,ni)
    set(gca,'FontSize',16)
    % hold on
    %  mesh(xedges,yedges,histmat')
    set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    xlabel('Px')
    ylabel('Py')
    
    xedges=[]; yedges=[]; histmat=[];
    xedges = linspace(-0.6e4,0.6e4,1e2); yedges =linspace(-0.6e4,0.6e4,1e2);
    histmat = hist2(eval(['xp_' int2str(ni)  '(aaa_' int2str(ni) ')']),eval(['yp_' int2str(ni)  '(aaa_' int2str(ni) ')']), xedges, yedges);
    
    subplot(2,3,ni+3)
    set(gca,'FontSize',16)
    % hold on
    %   mesh(xedges,yedges,histmat')
    set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    xlabel('Px')
    ylabel('Py')
    
    
end


%% transvers plot

l=10;

for ni=1:1:3
    
    eval(['x1_coor' int2str(ni) '=x_coor' int2str(ni) '+(xp_' int2str(ni) './zp_' int2str(ni) ').*l;']);
    eval(['y1_coor' int2str(ni) '=y_coor' int2str(ni) '+(yp_' int2str(ni) './zp_' int2str(ni) ').*l;']);
end

figure(ifig)
ifig=ifig+1;
suptitle(['after ' num2str(l) ' m; theta<' num2str(theta)...
    '; bandwith =' num2str(std(full_spectrum_1(aaa_1))/mean(full_spectrum_1(aaa_1)))])
for ni=1:1:3
    xedges=[]; yedges=[]; histmat=[];
    xedges = linspace(-400e-4,400e-4,1e2); yedges =linspace(-400e-4,400e-4,1e2);
    histmat = hist2(eval(['x1_coor' int2str(ni) ]),eval(['y1_coor' int2str(ni) ]), xedges, yedges);
    subplot(2,3,ni)
    set(gca,'FontSize',16)
%     mesh(xedges,yedges,histmat')
    set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    xlabel('x')
    ylabel('y')
    
    
    xedges=[]; yedges=[]; histmat=[];
    xedges = linspace(-400e-4,400e-4,1e2); yedges =linspace(-400e-4,400e-4,1e2);
    histmat = hist2(eval(['x1_coor' int2str(ni)  '(aaa_' int2str(ni) ')']),eval(['y1_coor' int2str(ni)  '(aaa_' int2str(ni) ')']), xedges, yedges);
    subplot(2,3,ni+3)
    set(gca,'FontSize',16)
%     mesh(xedges,yedges,histmat')
    set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    xlabel('x')
    ylabel('y')
    
end





stop
%%%% experimental part

figure(ifig)
ifig=ifig+1;
for ni=1:1:3
    subplot(1,3,ni)
    % scatter(eval(['x_coor' int2str(ni) ]),eval(['y_coor' int2str(ni) ]),10,eval(['Sx_' int2str(ni) ]), 'filled');
    
    % colorbar;
    % subplot(1,2,2)
    
    
    scatter(eval(['x_coor' int2str(ni) ]),...
        eval(['y_coor' int2str(ni) ]),5,...
        [eval(['Sx_' int2str(ni) ])...
        eval(['Sy_' int2str(ni) ])...
        eval(['Sz_' int2str(ni) ])], 'filled');
    colorbar;
    
    
    
    % scatter(eval(['x_coor' int2str(ni) '(aa_' int2str(ni) ')']),...
    %     eval(['y_coor' int2str(ni) '(aa_' int2str(ni) ')']),50,...
    %     [eval(['Sx_' int2str(ni) '(aa_' int2str(ni) ')'])...
    %      eval(['Sy_' int2str(ni) '(aa_' int2str(ni) ')'])...
    %      eval(['Sz_' int2str(ni) '(aa_' int2str(ni) ')'])], 'filled');
    %  colorbar;
    
    % set(gca,'FontSize',16)
    % xlabel('Sz')
    %  title({['std ' num2str(eval(['std(Sz_' int2str(ni)  '(aa_' int2str(ni) '))'])) ''],...
    %      ['mean ' num2str(eval(['mean(Sz_' int2str(ni)  '(aa_' int2str(ni) '))'])) '']})
    % suptitle(['after ' num2str(l) ' m'])
end
% filename = ['Stokes_' num2str(ifig) ];
% fname = [ filename '.png'];
% print('-dpng', fname);


figure(ifig)
ifig=ifig+1;
% image(x_coor1(aa_1),y_coor1(aa_1),Sx_1(aa_1));
% scatter(x_coor1(aa_1),y_coor1(aa_1), 10,Sx_1(aa_1), 'filled');
scatter(x_coor3(aa_3),y_coor3(aa_3), 10,[Sx_3(aa_3) Sy_3(aa_3) Sz_3(aa_3)]);
colorbar;






figure(ifig)
ifig=ifig+1;
subplot 211
hist(sqrt(Sx_1(aa_1).^2+Sy_1(aa_1).^2+Sz_1(aa_1).^2))
subplot 212
hist(sqrt(Sx_1.^2+Sy_1.^2+Sz_1.^2))

%%

% Phi=asin(xp_1./sqrt((xp_1.^2+yp_1.^2)));
Phi=atan(yp_1./xp_1);
% Snx_1=-Sx_1.*cos(2.*Phi)+Sz_1.*sin(2.*Phi);
Sny_1=Sy_1;
% Snz_1=-Sx_1.*sin(2.*Phi)-Sz_1.*cos(2.*Phi);
Snx_1=-Sx_1.*cos(2.*Phi)-Sz_1.*sin(2.*Phi);
Snz_1=Sx_1.*sin(2.*Phi)-Sz_1.*cos(2.*Phi);


%% Sx Sy Sz
figure(ifig)
ifig=ifig+1;
for ni=1:1:1
    subplot(1,3,1)
    hist(Snx_1, 10);
    xlabel('Sx')
    title({['std ' num2str(std(Snx_1)) ]; ['mean ' num2str(mean(Snx_1)) ]})
    subplot(1,3,2)
    hist(Sny_1, 10);
    xlabel('Sy')
    title({['std ' num2str(std(Sny_1)) ]; ['mean ' num2str(mean(Sny_1)) ]})
    % suptitle(['after ' num2str(l) ' m'])
    subplot(1,3,3)
    hist(Snz_1, 10);
    xlabel('Sz')
    title({['std ' num2str(std(Snz_1)) ]; ['mean ' num2str(mean(Snz_1)) ]})
end

figure(ifig)
ifig=ifig+1;
% hist((xp_1.^2+yp_1.^2)./(yp_1.^2),100)
hist(Phi)

% Phi(aa_1)=asin(xp_1(aa_1)./sqrt((xp_1(aa_1).^2+yp_1(aa_1).^2)));
%
% Snx_1(aa_1)=-Sx_1(aa_1).*cos(2.*Phi(aa_1))+Sz_1(aa_1).*sin(2.*Phi(aa_1));
% Sny_1(aa_1)=Sy_1(aa_1);
% Snz_1(aa_1)=-Sx_1(aa_1).*sin(2.*Phi(aa_1))+Sz_1(aa_1).*cos(2.*Phi(aa_1));

%% Sx Sy Sz
figure(ifig)
ifig=ifig+1;
for ni=1:1:1
    subplot(1,3,1)
    hist(Snx_1(aa_1), 10);
    title({['std ' num2str(std(Snx_1(aa_1))) ]; ['mean ' num2str(mean(Snx_1(aa_1))) ]})
    xlabel('Sx')
    subplot(1,3,2)
    hist(Sny_1(aa_1), 10);
    xlabel('Sy')
    title({['std ' num2str(std(Sny_1(aa_1))) ]; ['mean ' num2str(mean(Sny_1(aa_1))) ]})
    subplot(1,3,3)
    hist(Snz_1(aa_1), 10);
    xlabel('Sz')
    title({['std ' num2str(std(Snz_1(aa_1))) ]; ['mean ' num2str(mean(Snz_1(aa_1))) ]})
end

figure(ifig)
ifig=ifig+1;
subplot 131
plot(phot_angle_1,Sx_1,'.r')
subplot 132
plot(phot_angle_1,Sy_1,'.r')
subplot 133
plot(phot_angle_1,Sz_1,'.r')
% hist(sqrt((xp_1.^2+yp_1.^2))./yp_1,100)




figure(ifig)
ifig=ifig+1;
subplot 131
plot(phot_angle_1(aa_1),Sx_1(aa_1),'.r')
subplot 132
plot(phot_angle_1(aa_1),Sy_1(aa_1),'.r')
subplot 133
plot(phot_angle_1(aa_1),Sz_1(aa_1),'.r')













