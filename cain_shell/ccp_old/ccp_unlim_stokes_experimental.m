% one beam several events
clear all; close all; clc;

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


%STOKES(1,:)=[1 0 0];
%STOKES(2,:)=[0 1 0];
%STOKES(3,:)=[0 0 1];

aaa=[250 275 225];
for var_for_scan=aaa;
    qq=qq+1;
    
 
    %rflags.STOKES=STOKES(var_for_scan,:);
    rflags.pulseE=0.2; %laser puse energy [J]
    
    beam_phasespace=dlmread(['Anna_beams/eli_lowen_oned_wp_20_11_run00' num2str(qq) '.w5.sdds.asci'],'',10,0);;
    beam_phasespace(:,6)=beam_phasespace(:,6)-114.7/0.511;
    defoc_param=0.9;%defocusing parameter by defoult 1 is no defocusing
    add_param=std(beam_phasespace(:,6));
BASE_DIRECTORY = [pwd '/Echarge_' num2str(var_for_scan) '_pulseE_' num2str(rflags.pulseE) '/'];
    %%
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
filename = ['EXP_plot_' num2str(add_param) '_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);






theta=1.2e-4;
ww=0;
for ni=aaa
    ww=ww+1
    eval(['aa_' num2str(ww) '=find(abs(phot_angle_' num2str(ww) ')<' num2str(theta) ');']);
    num_in_b(ww)=eval(['length(full_spectrum_' int2str(ww) '(aa_' int2str(ww) '));']);
    
end

nbin=50;
min_val_linspase=100
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
filename = ['EXP_plot_' num2str(add_param) '_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);


%% Stokes 

figure(ifig)
ifig=ifig+1;
for ni=1:1:3


subplot(3,3,ni)
hist(eval(['Sx_' int2str(ni)  '(aa_' int2str(ni) ')']), 50);
set(gca,'FontSize',16)
xlabel('Sx')
 title({['std ' num2str(eval(['std(Sx_' int2str(ni)  '(aa_' int2str(ni) '))'])) ''],...
     ['mean ' num2str(eval(['mean(Sx_' int2str(ni)  '(aa_' int2str(ni) '))'])) '']})

subplot(3,3,ni+3)
hist(eval(['Sy_' int2str(ni) '(aa_' int2str(ni) ')']), 100);
set(gca,'FontSize',16)
xlabel('Sy')
 title({['std ' num2str(eval(['std(Sy_' int2str(ni)  '(aa_' int2str(ni) '))'])) ''],...
     ['mean ' num2str(eval(['mean(Sy_' int2str(ni)  '(aa_' int2str(ni) '))'])) '']})
% suptitle(['after ' num2str(l) ' m'])

subplot(3,3,ni+6)
hist(eval(['Sz_' int2str(ni) '(aa_' int2str(ni) ')']), 100);
set(gca,'FontSize',16)
xlabel('Sz')
 title({['std ' num2str(eval(['std(Sz_' int2str(ni)  '(aa_' int2str(ni) '))'])) ''],...
     ['mean ' num2str(eval(['mean(Sz_' int2str(ni)  '(aa_' int2str(ni) '))'])) '']})
 
end
suptitle(['in theta<' num2str(theta) ''])
filename = ['EXP_plot_' num2str(add_param) '_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);


figure(ifig)
ifig=ifig+1;
for ni=1:1:3

% set(gca,'FontSize',16)
subplot(3,3,ni)
hist(eval(['Sx_' int2str(ni)]), 50);
% set(gca,'FontSize',16)
xlabel('Sx')
 title({['std ' num2str(eval(['std(Sx_' int2str(ni) ')'])) ''],...
     ['mean ' num2str(eval(['mean(Sx_' int2str(ni) ')'])) '']})
%  title({['Sx; std ' num2str(eval(['std(Sx_' int2str(ni) ')'])) '; mean ' num2str(eval(['mean(Sx_' int2str(ni) ')'])) '']})

subplot(3,3,ni+3)
hist(eval(['Sy_' int2str(ni)]), 100);
% set(gca,'FontSize',16)
xlabel('Sy')
 title({['std ' num2str(eval(['std(Sy_' int2str(ni) ')'])) ''],...
     ['mean ' num2str(eval(['mean(Sy_' int2str(ni) ')'])) '']})
%  title({['Sy; std ' num2str(eval(['std(Sy_' int2str(ni) ')'])) '; mean ' num2str(eval(['mean(Sy_' int2str(ni) ')'])) '']})
 
subplot(3,3,ni+6)
hist(eval(['Sz_' int2str(ni)]), 100);
% set(gca,'FontSize',16)
 xlabel('Sz')
title({['std ' num2str(eval(['std(Sz_' int2str(ni) ')'])) ''],...
     [' mean ' num2str(eval(['mean(Sz_' int2str(ni) ')'])) '']})
 
end
suptitle(['all photons'])
filename = ['12121EXP_plot_' num2str(add_param) '_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng','-r400', fname);
%%%%%%%%%%%%%%%%%


figure(ifig)
ifig=ifig+1;
for ni=1:1:3


subplot(1,3,ni)
hist(eval(['Sx_' int2str(ni)  '(aa_' int2str(ni) ')']), 50);
%set(gca,'FontSize',16)
xlabel('Sx')
 title({['std ' num2str(eval(['std(Sx_' int2str(ni)  '(aa_' int2str(ni) '))'])) ''],...
     ['mean ' num2str(eval(['mean(Sx_' int2str(ni)  '(aa_' int2str(ni) '))'])) '']})
end
filename = ['Stokes_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);



figure(ifig)
ifig=ifig+1;
for ni=1:1:3


subplot(1,3,ni)
hist(eval(['Sy_' int2str(ni) '(aa_' int2str(ni) ')']), 100);
%set(gca,'FontSize',16)
xlabel('Sy')
 title({['std ' num2str(eval(['std(Sy_' int2str(ni)  '(aa_' int2str(ni) '))'])) ''],...
     ['mean ' num2str(eval(['mean(Sy_' int2str(ni)  '(aa_' int2str(ni) '))'])) '']})
% suptitle(['after ' num2str(l) ' m'])
end
filename = ['Stokes_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);


figure(ifig)
ifig=ifig+1;
for ni=1:1:3


subplot(1,3,ni)
hist(eval(['Sz_' int2str(ni) '(aa_' int2str(ni) ')']), 100);
%set(gca,'FontSize',16)
xlabel('Sz')
 title({['std ' num2str(eval(['std(Sz_' int2str(ni)  '(aa_' int2str(ni) '))'])) ''],...
     ['mean ' num2str(eval(['mean(Sz_' int2str(ni)  '(aa_' int2str(ni) '))'])) '']})
% suptitle(['after ' num2str(l) ' m'])
end
filename = ['Stokes_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);

stop

%% Stokes  temperature plot
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
suptitle(['all photons'])
filename = ['EXP_plot_' num2str(add_param) '_' num2str(ifig) ];
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
filename = ['EXP_plot_' num2str(add_param) '_' num2str(ifig) ];
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
%     mesh(xedges,yedges,histmat')
     set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    xlabel('Px')
    ylabel('Py')
    
    xedges=[]; yedges=[]; histmat=[];
    xedges = linspace(-0.6e4,0.6e4,1e2); yedges =linspace(-0.6e4,0.6e4,1e2);
    histmat = hist2(eval(['xp_' int2str(ni)  '(aaa_' int2str(ni) ')']),eval(['yp_' int2str(ni)  '(aaa_' int2str(ni) ')']), xedges, yedges);
    
    subplot(2,3,ni+3)
    set(gca,'FontSize',16)
    % hold on
%        mesh(xedges,yedges,histmat')
    set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    xlabel('Px')
    ylabel('Py')
    
    
end
filename = ['EXP_plot_' num2str(add_param) '_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);

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
    mesh(xedges,yedges,histmat')
    % set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    xlabel('x')
    ylabel('y')
    
    
    xedges=[]; yedges=[]; histmat=[];
    xedges = linspace(-400e-4,400e-4,1e2); yedges =linspace(-400e-4,400e-4,1e2);
    histmat = hist2(eval(['x1_coor' int2str(ni)  '(aaa_' int2str(ni) ')']),eval(['y1_coor' int2str(ni)  '(aaa_' int2str(ni) ')']), xedges, yedges);
    subplot(2,3,ni+3)
    set(gca,'FontSize',16)
    mesh(xedges,yedges,histmat')
    % set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    xlabel('x')
    ylabel('y')
    
end
filename = ['EXP_plot_' num2str(add_param) '_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);


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
    %mesh(xedges,yedges,histmat')
    set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    xlabel('x')
    ylabel('y')
    
    
    xedges=[]; yedges=[]; histmat=[];
    xedges = linspace(-10e-4,10e-4,1e2); yedges =linspace(-10e-4,10e-4,1e2);
    histmat = hist2(eval(['x1_coor' int2str(ni)  '(aaa_' int2str(ni) ')']),eval(['y1_coor' int2str(ni)  '(aaa_' int2str(ni) ')']), xedges, yedges);
    subplot(2,3,ni+3)
    set(gca,'FontSize',16)
    %mesh(xedges,yedges,histmat')
    set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    xlabel('x')
    ylabel('y')
    
end
filename = ['EXP_plot_' num2str(add_param) '_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);

%% 

%% after propagation


bu_factor=5
theta=bu_factor*5e-5;
num_s =3;
w_size=bu_factor*5e-5;
nx_bin=50;
ny_bin=50;
xcor_allb=[]; ycor_allb=[];
    eval(['x1_coor' int2str(num_s) '=[];']);
    eval(['y1_coor' int2str(num_s) '=[];']);
eval(['aa_' num2str(num_s) '=find(abs(phot_angle_' num2str(num_s) ')<' num2str(theta) ');']);

for l=1%[0:0.0001:0.01]% 0.1:0.1:0.5];%
   
    
    eval(['ycor_nop=y_coor' int2str(num_s) '(aa_' int2str(num_s) ');']);
    eval(['xcor_nop=x_coor' int2str(num_s) '(aa_' int2str(num_s) ');']);
    eval(['Sz=Sz_' int2str(num_s) ';']);
    eval(['Sz_sep=Sz_' int2str(num_s) '(aa_' int2str(num_s) ');']);
    eval(['x1_coor' int2str(num_s) '=x_coor' int2str(num_s) '+(xp_' int2str(num_s) './zp_' int2str(num_s) ').*l;']);
    eval(['y1_coor' int2str(num_s) '=y_coor' int2str(num_s) '+(yp_' int2str(num_s) './zp_' int2str(num_s) ').*l;']);
    eval(['xcor_allb=x1_coor' int2str(num_s) ';']);
    eval(['ycor_allb=y1_coor' int2str(num_s) ';']);
    eval(['xcor_n=xcor_nop+(xp_' int2str(num_s) '(aa_' int2str(num_s) ')./zp_' int2str(num_s) '(aa_' int2str(num_s) ')).*l;']);
    eval(['ycor_n=ycor_nop+(yp_' int2str(num_s) '(aa_' int2str(num_s) ')./zp_' int2str(num_s) '(aa_' int2str(num_s) ')).*l;']);
    
    d1=zeros(ny_bin,nx_bin);
    d2=zeros(ny_bin,nx_bin);
    f=zeros(ny_bin,nx_bin);
    d1_theta=zeros(ny_bin,nx_bin);
    d2_theta=zeros(ny_bin,nx_bin);
    f_theta=zeros(ny_bin,nx_bin);
    
    x_l=2*w_size/nx_bin;
    y_l=2*w_size/ny_bin;
    
    
    for ni=1:1:length(xcor_allb);
        if (xcor_allb(ni)<w_size && xcor_allb(ni)>-w_size)
            if (ycor_allb(ni)<w_size && ycor_allb(ni)>-w_size)
                d1(1+floor((ycor_allb(ni)+w_size)./y_l),1+floor((xcor_allb(ni)+w_size)./x_l))=...
                    d1(1+floor((ycor_allb(ni)+w_size)./y_l),1+floor((xcor_allb(ni)+w_size)./x_l))+Sz(ni);
                d2(1+floor((ycor_allb(ni)+w_size)./y_l),1+floor((xcor_allb(ni)+w_size)./x_l))=...
                    d2(1+floor((ycor_allb(ni)+w_size)./y_l),1+floor((xcor_allb(ni)+w_size)./x_l))+1;
            end
        end
    end
    
    
    for nyi=1:1:ny_bin
        for nxi=1:1:nx_bin
            if(d2(nyi,nxi)==0)
                f(nyi,nxi)=0;
            else
                f(nyi,nxi)=d1(nyi,nxi)./d2(nyi,nxi);
            end
        end
    end
    
    
    
    for ni=1:1:length(xcor_n);
        if (xcor_n(ni)<w_size && xcor_n(ni)>-w_size)
            if (ycor_n(ni)<w_size && ycor_n(ni)>-w_size)
                d1_theta(1+floor((ycor_n(ni)+w_size)./y_l),1+floor((xcor_n(ni)+w_size)./x_l))=d1_theta(1+floor((ycor_n(ni)+w_size)./y_l),1+floor((xcor_n(ni)+w_size)./x_l))+Sz_sep(ni);
                d2_theta(1+floor((ycor_n(ni)+w_size)./y_l),1+floor((xcor_n(ni)+w_size)./x_l))=d2_theta(1+floor((ycor_n(ni)+w_size)./y_l),1+floor((xcor_n(ni)+w_size)./x_l))+1;
            end
        end
    end
    for nyi=1:1:ny_bin
        for nxi=1:1:nx_bin
            if(d2_theta(nyi,nxi)==0)
                f_theta(nyi,nxi)=0;
            else
                f_theta(nyi,nxi)=d1_theta(nyi,nxi)./d2_theta(nyi,nxi);
            end
        end
    end
    
    
    
    
    
    
    
    
    figure(100)
    % ifig=ifig+1;
    subplot 221
    set(pcolor(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        f),'EdgeColor','none')
    colorbar('SouthOutside')
    title('prop')
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    
    subplot 222
    set(pcolor(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        d2),'EdgeColor','none')
    colorbar('SouthOutside')
    title('prop')
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    suptitle(['L=' num2str(l)])
    
    
    subplot 223
    set(pcolor(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        f_theta),'EdgeColor','none')
    colorbar('SouthOutside')
    title('prop')
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    
    subplot 224
    set(pcolor(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        d2_theta),'EdgeColor','none')
    colorbar('SouthOutside')
    title('prop')
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    suptitle(['L=' num2str(l) ' m; Stokes=' num2str(STOKES(num_s,:)) ])
    
    
    filename = ['EXP_plot_' num2str(1) '_' num2str(l*1e6) ];
    fname = [ filename '.png'];
    print('-dpng', fname);
    
end












