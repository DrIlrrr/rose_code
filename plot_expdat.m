% one beam several events

clear all; close all; clc;


old=load(['old.dat']);
new=load(['new.dat']);

%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss
% WEIGHT=3;
TIME_COORDINATE=4;%now we start use s(m) for caine
X_COORDINATE=5;
Y_COORDINATE=6;
ENERGY_OF_PARTICLE=8;
X_MOMENTUM=9;
Y_MOMENTUM=10;
S_MOMENTUM=11;
% % POLARISATION: 12 13 14
% N_COMPTON_HIT=15;
% TURN_LAST_COMPTON_HIT=16;
ifig=0;

NumberMP_old = length(old);
NumberMP_new = length(new);
Echarge = 1.60e-19;% Charge of electron [c]
Ne = 250e-12/Echarge;% Number electrons in bunch
weight_old = Ne/NumberMP_old;
weight_new = Ne/NumberMP_new;




nbin_plot=50;

nbin_old=floor((max(old(:,5))+min(old(:,5)))*1e7);% nbin to make spectrum per m
nbin_new=floor((max(new(:,5))+min(new(:,5)))*1e7);% nbin to make spectrum per m
xs_old=linspace(min(old(:,5)),max(old(:,5)),nbin_old);
xs_new=linspace(min(new(:,5)),max(new(:,5)),nbin_new);


ifig=ifig+1;
figure(ifig)
subplot 211
bar(xs_old,hist(old(:,5),nbin_old)*weight_old,'grouped','hist','g')
subplot 212
bar(xs_new,hist(new(:,5),nbin_new)*weight_new,'grouped','hist','b')
filename = ['comp_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);


ifig=ifig+1;
figure(ifig)
hold on
bar(xs_old,hist(old(:,5),nbin_old)*weight_old,'grouped','hist','g')
bar(xs_new,hist(new(:,5),nbin_new)*weight_new,'grouped','hist','b')
hold off
filename = ['comp_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);



ifig=ifig+1;
figure(ifig)
hold on
plot(xs_old,hist(old(:,5),nbin_old)*weight_old,'g')
plot(xs_new,hist(new(:,5),nbin_new)*weight_new,'b')
hold off
filename = ['comp_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);



%% XY z at IP
    
    xedges = linspace(-0.6e-4,0.6e-4,0.5e2); yedges =linspace(-0.6e-4,0.6e-4,0.5e2);
    histmat_old = hist2(old(:,5),old(:,6), xedges, yedges);
    histmat_new = hist2(new(:,5),new(:,6), xedges, yedges);
    
    figure(ifig)
    ifig=ifig+1;
    
    subplot 121
    % mesh(xedges,yedges,histmat')
    set(pcolor(xedges,yedges,histmat_old'),'EdgeColor','none')
    colormap(jet)
 
%      title({['old '];['std x=' num2str(std(old(:,5))*1e6) ' [\mu m]; mean of x=' num2str(mean(old(:,5))*1e6) ' [\mu m]'];...
%         ['std y=' num2str(std(old(:,6))*1e6) ' [\mu m]; mean of y=' num2str(mean(old(:,6))*1e6)  '[\mu m]']});
set(gca,'FontSize',12)  
title({['old '];['std x=' num2str(std(old(:,5))*1e6) ' [\mu m]; '];...
    ['std y=' num2str(std(old(:,6))*1e6) ' [\mu m];'];[]});
   set(gca,'FontSize',16)
     xlabel('x')
    ylabel('y')
     subplot 122
 % mesh(xedges,yedges,histmat')
    set(pcolor(xedges,yedges,histmat_new'),'EdgeColor','none')
    colormap(jet)
set(gca,'FontSize',12)
    title({['new '];['std x=' num2str(std(new(:,5))*1e6) ' [\mu m]; '];...
        ['std y=' num2str(std(new(:,6))*1e6) ' [\mu m];'];[]});
set(gca,'FontSize',16)
    xlabel('x')
    ylabel('y')    
    
    filename = ['comp_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);






ifig=ifig+1;
figure(ifig)
subplot 211
bar(hist(old(:,6),nbin_old)*weight_old,'grouped','hist','g')
subplot 212
bar(hist(new(:,6),nbin_new)*weight_new,'grouped','hist','b')
filename = ['comp_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);







%% XY z at IP
    
    xedges = linspace(-0.6e-4,0.6e-4,0.5e2); yedges =linspace(-0.6e-4,0.6e-4,0.5e2);
    histmat_old = hist2(old(:,5),old(:,6), xedges, yedges);
    histmat_new = hist2(new(:,5),new(:,6), xedges, yedges);
    
    figure(ifig)
    ifig=ifig+1;
    
    subplot 121
     mesh(xedges,yedges,histmat_old')
%     set(pcolor(xedges,yedges,histmat_old'),'EdgeColor','none')
    colormap(jet)
 
%      title({['old '];['std x=' num2str(std(old(:,5))*1e6) ' [\mu m]; mean of x=' num2str(mean(old(:,5))*1e6) ' [\mu m]'];...
%         ['std y=' num2str(std(old(:,6))*1e6) ' [\mu m]; mean of y=' num2str(mean(old(:,6))*1e6)  '[\mu m]']});
set(gca,'FontSize',12)  
title({['old '];['std x=' num2str(std(old(:,5))*1e6) ' [\mu m]; '];...
    ['std y=' num2str(std(old(:,6))*1e6) ' [\mu m];'];[]});
   set(gca,'FontSize',16)
     xlabel('x')
    ylabel('y')
     subplot 122
  mesh(xedges,yedges,histmat_new')
%     set(pcolor(xedges,yedges,histmat_new'),'EdgeColor','none')
    colormap(jet)
set(gca,'FontSize',12)
    title({['new '];['std x=' num2str(std(new(:,5))*1e6) ' [\mu m]; '];...
        ['std y=' num2str(std(new(:,6))*1e6) ' [\mu m];'];[]});
set(gca,'FontSize',16)
    xlabel('x')
    ylabel('y')    
    
    filename = ['comp_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);





%% XY z at IP
    aa=0.04e6;
    xedges = linspace(-aa,aa,0.5e2); yedges =linspace(-aa,aa,0.5e2);
    histmat_old = hist2(old(:,9),old(:,10), xedges, yedges);
    histmat_new = hist2(new(:,9),new(:,10), xedges, yedges);
    
    figure(ifig)
    ifig=ifig+1;
    
    subplot 121
%      mesh(xedges,yedges,histmat_old')
    set(pcolor(xedges,yedges,histmat_old'),'EdgeColor','none')
    colormap(jet)
 
%      title({['old '];['std x=' num2str(std(old(:,5))*1e6) ' [\mu m]; mean of x=' num2str(mean(old(:,5))*1e6) ' [\mu m]'];...
%         ['std y=' num2str(std(old(:,6))*1e6) ' [\mu m]; mean of y=' num2str(mean(old(:,6))*1e6)  '[\mu m]']});
set(gca,'FontSize',12)  
title({['old '];['std x=' num2str(std(old(:,5))*1e6) ' [\mu m]; '];...
    ['std y=' num2str(std(old(:,6))*1e6) ' [\mu m];'];[]});
   set(gca,'FontSize',16)
     xlabel('P_x')
    ylabel('P_y')
     subplot 122
%   mesh(xedges,yedges,histmat_new')
   set(pcolor(xedges,yedges,histmat_new'),'EdgeColor','none')
    colormap(jet)
set(gca,'FontSize',12)
    title({['new '];['std x=' num2str(std(new(:,5))*1e6) ' [\mu m]; '];...
        ['std y=' num2str(std(new(:,6))*1e6) ' [\mu m];'];[]});
set(gca,'FontSize',16)
    xlabel('P_x')
    ylabel('P_y')    
    
    filename = ['comp_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);




%% XY z at IP
    aa=0.04e6;
    xedges = linspace(-aa,aa,0.5e2); yedges =linspace(-aa,aa,0.5e2);
    histmat_old = hist2(old(:,9),old(:,10), xedges, yedges);
    histmat_new = hist2(new(:,9),new(:,10), xedges, yedges);
    
    figure(ifig)
    ifig=ifig+1;
    
    subplot 121
      mesh(xedges,yedges,histmat_old')
%     set(pcolor(xedges,yedges,histmat_old'),'EdgeColor','none')
    colormap(jet)
 
%      title({['old '];['std x=' num2str(std(old(:,5))*1e6) ' [\mu m]; mean of x=' num2str(mean(old(:,5))*1e6) ' [\mu m]'];...
%         ['std y=' num2str(std(old(:,6))*1e6) ' [\mu m]; mean of y=' num2str(mean(old(:,6))*1e6)  '[\mu m]']});
set(gca,'FontSize',12)  
title({['old '];['std x=' num2str(std(old(:,9))) ' [\mu m]; '];...
    ['std y=' num2str(std(old(:,10))) ' [\mu m];'];[]});
   set(gca,'FontSize',16)
     xlabel('P_x')
    ylabel('P_y')
     subplot 122
   mesh(xedges,yedges,histmat_new')
%    set(pcolor(xedges,yedges,histmat_new'),'EdgeColor','none')
    colormap(jet)
set(gca,'FontSize',12)
    title({['new '];['std x=' num2str(std(new(:,9))) ' [\mu m]; '];...
        ['std y=' num2str(std(new(:,10))) ' [\mu m];'];[]});
set(gca,'FontSize',16)
    xlabel('P_x')
    ylabel('P_y')    
    
    filename = ['comp_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);


%% XY z at IP
    aa=0.04e6;
    xedges = linspace(-aa,aa,0.5e2); yedges =linspace(-0.6e-4,0.6e-4,0.5e2);
    histmat_old = hist2(old(:,10),old(:,6), xedges, yedges);
    histmat_new = hist2(new(:,10),new(:,6), xedges, yedges);
    
    figure(ifig)
    ifig=ifig+1;
    
    subplot 121
%      mesh(xedges,yedges,histmat_old')
    set(pcolor(xedges,yedges,histmat_old'),'EdgeColor','none')
    colormap(jet)
 
%      title({['old '];['std x=' num2str(std(old(:,5))*1e6) ' [\mu m]; mean of x=' num2str(mean(old(:,5))*1e6) ' [\mu m]'];...
%         ['std y=' num2str(std(old(:,6))*1e6) ' [\mu m]; mean of y=' num2str(mean(old(:,6))*1e6)  '[\mu m]']});
set(gca,'FontSize',12)  
title({['old '];['std x=' num2str(std(old(:,5))*1e6) ' [\mu m]; '];...
    ['std y=' num2str(std(old(:,6))*1e6) ' [\mu m];'];[]});
   set(gca,'FontSize',16)
     xlabel('y')
    ylabel('P_y')
     subplot 122
%   mesh(xedges,yedges,histmat_new')
   set(pcolor(xedges,yedges,histmat_new'),'EdgeColor','none')
    colormap(jet)
set(gca,'FontSize',12)
    title({['new '];['std x=' num2str(std(new(:,5))*1e6) ' [\mu m]; '];...
        ['std y=' num2str(std(new(:,6))*1e6) ' [\mu m];'];[]});
set(gca,'FontSize',16)
    xlabel('y')
    ylabel('P_y')    
    
    filename = ['comp_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);


%% XY z at IP
    aa=0.04e6;
    xedges = linspace(-aa,aa,0.5e2); yedges =linspace(-0.6e-4,0.6e-4,0.5e2);
    histmat_old = hist2(old(:,9),old(:,5), xedges, yedges);
    histmat_new = hist2(new(:,9),new(:,5), xedges, yedges);
    
    figure(ifig)
    ifig=ifig+1;
    
    subplot 121
%      mesh(xedges,yedges,histmat_old')
    set(pcolor(xedges,yedges,histmat_old'),'EdgeColor','none')
    colormap(jet)
 
%      title({['old '];['std x=' num2str(std(old(:,5))*1e6) ' [\mu m]; mean of x=' num2str(mean(old(:,5))*1e6) ' [\mu m]'];...
%         ['std y=' num2str(std(old(:,6))*1e6) ' [\mu m]; mean of y=' num2str(mean(old(:,6))*1e6)  '[\mu m]']});
set(gca,'FontSize',12)  
title({['old '];['std x=' num2str(std(old(:,5))*1e6) ' [\mu m]; '];...
    ['std y=' num2str(std(old(:,6))*1e6) ' [\mu m];'];[]});
   set(gca,'FontSize',16)
     xlabel('x')
    ylabel('P_x')
     subplot 122
%   mesh(xedges,yedges,histmat_new')
   set(pcolor(xedges,yedges,histmat_new'),'EdgeColor','none')
    colormap(jet)
set(gca,'FontSize',12)
    title({['new '];['std x=' num2str(std(new(:,5))*1e6) ' [\mu m]; '];...
        ['std y=' num2str(std(new(:,6))*1e6) ' [\mu m];'];[]});
set(gca,'FontSize',16)
    xlabel('x')
    ylabel('P_x')    
    
    filename = ['comp_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);

nbin_old=floor((max(old(:,X_MOMENTUM))-min(old(:,X_MOMENTUM)))/1e3);% nbin to make spectrum per m
nbin_new=floor((max(new(:,X_MOMENTUM))-min(new(:,X_MOMENTUM)))/1e3);% nbin to make spectrum per m
xs_old=linspace(min(old(:,X_MOMENTUM)),max(old(:,X_MOMENTUM)),nbin_old);
xs_new=linspace(min(new(:,X_MOMENTUM)),max(new(:,X_MOMENTUM)),nbin_new);



ifig=ifig+1;
figure(ifig)
hold on
bar(xs_old,hist(old(:,X_MOMENTUM),nbin_old)*weight_old,'grouped','hist','g')
bar(xs_new,hist(new(:,X_MOMENTUM),nbin_new)*weight_new,'grouped','hist','b')
hold off
xlabel('P_x')
filename = ['comp_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);



ifig=ifig+1;
figure(ifig)
hold on
plot(xs_old,hist(old(:,X_MOMENTUM),nbin_old)*weight_old,'g')
plot(xs_new,hist(new(:,X_MOMENTUM),nbin_new)*weight_new,'b')
hold off
xlabel('P_x')
legend('old','new',0)
filename = ['comp_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);


%%

nbin_old=floor((max(old(:,Y_MOMENTUM))-min(old(:,Y_MOMENTUM)))/1e3);% nbin to make spectrum per m
nbin_new=floor((max(new(:,Y_MOMENTUM))-min(new(:,Y_MOMENTUM)))/1e3);% nbin to make spectrum per m
xs_old=linspace(min(old(:,Y_MOMENTUM)),max(old(:,Y_MOMENTUM)),nbin_old);
xs_new=linspace(min(new(:,Y_MOMENTUM)),max(new(:,Y_MOMENTUM)),nbin_new);



ifig=ifig+1;
figure(ifig)
hold on
bar(xs_old,hist(old(:,Y_MOMENTUM),nbin_old)*weight_old,'grouped','hist','g')
bar(xs_new,hist(new(:,Y_MOMENTUM),nbin_new)*weight_new,'grouped','hist','b')
hold off
xlabel('P_y')
filename = ['comp_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);



ifig=ifig+1;
figure(ifig)
hold on
plot(xs_old,hist(old(:,Y_MOMENTUM),nbin_old)*weight_old,'g')
plot(xs_new,hist(new(:,Y_MOMENTUM),nbin_new)*weight_new,'b')
hold off
xlabel('P_y')
legend('old','new',0)
filename = ['comp_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);


%%

nbin_old=floor((max(old(:,ENERGY_OF_PARTICLE))-min(old(:,ENERGY_OF_PARTICLE)))/1e4);% nbin to make spectrum per m
nbin_new=floor((max(new(:,ENERGY_OF_PARTICLE))-min(new(:,ENERGY_OF_PARTICLE)))/1e4);% nbin to make spectrum per m
xs_old=linspace(min(old(:,ENERGY_OF_PARTICLE)),max(old(:,ENERGY_OF_PARTICLE)),nbin_old);
xs_new=linspace(min(new(:,ENERGY_OF_PARTICLE)),max(new(:,ENERGY_OF_PARTICLE)),nbin_new);



ifig=ifig+1;
figure(ifig)
hold on
bar(xs_old,hist(old(:,ENERGY_OF_PARTICLE),nbin_old)*weight_old,'grouped','hist','g')
bar(xs_new,hist(new(:,ENERGY_OF_PARTICLE),nbin_new)*weight_new,'grouped','hist','b')
hold off
xlabel('E')
filename = ['comp_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);



ifig=ifig+1;
figure(ifig)
hold on
plot(xs_old,hist(old(:,ENERGY_OF_PARTICLE),nbin_old)*weight_old,'g')
plot(xs_new,hist(new(:,ENERGY_OF_PARTICLE),nbin_new)*weight_new,'b')
hold off
xlabel('E')
legend('old','new',0)
filename = ['comp_fig_' num2str(ifig)];
fname = [ filename '.png'];
print('-dpng', fname);

%%



pxold=old(:,Y_MOMENTUM);
pxnew=new(:,Y_MOMENTUM);


av_old=mean(pxold.^4)
av_new=mean(pxnew.^4)
rms_old=std(pxold.^4)
rms_new=std(pxnew.^4)

av_old=mean(pxold.^2)
av_new=mean(pxnew.^2)
rms_old=std(pxold.^2)
rms_new=std(pxnew.^2)


