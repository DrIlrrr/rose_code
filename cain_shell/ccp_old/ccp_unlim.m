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


aaa=[250 275 225 2500];
rflags.pulseE=0.2;

for var_for_scan=aaa
    qq=qq+1;
    
    defoc_param=1;%defocusing parameter by defoult 1 is no defocusing
    BASE_DIRECTORY = [pwd '/Echarge_' num2str(var_for_scan) '_pulseE_' num2str(rflags.pulseE) '/'];
        
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
    eval(['weigth_' int2str(qq) '=weigth;'])
    
end


color_line={'-r','-b','-g','-m','-y'};
color_dot={'-r','-b','-g','-m','-y'};






figure(ifig)
ifig=ifig+1;
for ni=1:1:qq
hold on
plot(linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))),...
    histc(eval(['full_spectrum_' int2str(ni) ]),linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))))*eval(['weigth_' int2str(ni) ]),...
    color_line{ni},'LineWidth',2)

% plot(linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))),...
%     histc(eval(['full_spectrum_' int2str(ni) ]),linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))))*weigth,...
%     'LineWidth',2);

grid on
hold off
end
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')
legend('250nC','275nC','225nC',0)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);


ww=0;
for ni=aaa
    ww=ww+1
    eval(['aa_' num2str(ww) '=find(abs(phot_angle_' num2str(ww) '));']);
    total_num_in_b(ww)=eval(['length(full_spectrum_' int2str(ww) '(aa_' int2str(ww) '))*weigth_' int2str(ww) ';']);
    
end

theta=18.2e-5;
ww=0;
for ni=aaa
    ww=ww+1
    eval(['aa_' num2str(ww) '=find(abs(phot_angle_' num2str(ww) ')<' num2str(theta) ');']);
    num_in_b(ww)=eval(['length(full_spectrum_' int2str(ww) '(aa_' int2str(ww) '))*weigth_' int2str(ww) ';']);
    
end


figure(ifig)
ifig=ifig+1;
subplot 211
hold on
plot(aaa,total_num_in_b,'-.xr','LineWidth',0.5)
plot(aaa,total_num_in_b,'xr','LineWidth',3)
hold off
grid on
title('total')
ylim([0 max(total_num_in_b)+max(total_num_in_b)*1e-1])
set(gca,'FontSize',16)
ylabel('number of scattered photons')
subplot 212
hold on
plot(aaa,num_in_b,'-.xb','LineWidth',0.5)
plot(aaa,num_in_b,'xb','LineWidth',3)
title(['in Theta=' num2str(theta) ])
hold off
grid on
ylim([0 max(num_in_b)+max(num_in_b)*1e-1])
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('Charge')
% legend('\delta t=0 [ps]','\delta t=1 [ps]','\delta t=2 [ps]',0)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);





nbin=50;
min_val_linspase=2500
figure(ifig)
ifig=ifig+1;
hold on
cc=jet(ww);
for ni=1:1:ww
plot(linspace(min_val_linspase,max(eval(['full_spectrum_' int2str(ni) ])),nbin),...
    histc(eval(['full_spectrum_' int2str(ni) '(aa_' int2str(ni) ')']),...
    linspace(min_val_linspase,max(eval(['full_spectrum_' int2str(ni) ])),nbin))*eval(['weigth_' int2str(ni) ]),'color',cc(ni,:),'LineWidth',2)
end
grid on
set(gca,'FontSize',16)
title({['Spectrum of scattered photons in Theta=' num2str(theta) ' [rad]'],...
     ['bandwith =' num2str(std(full_spectrum_1(aa_1))/mean(full_spectrum_1(aa_1))) ...
     '; ' num2str(std(full_spectrum_2(aa_2))/mean(full_spectrum_2(aa_2))) ';'...
      num2str(std(full_spectrum_3(aa_3))/mean(full_spectrum_3(aa_3))) '']})
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')
legend(['250nC; photons = ' num2str(num_in_b(1),'%2.1e') ],...
    ['275nC; photons = ' num2str(num_in_b(2),'%2.1e') ],...
    ['225nC; photons = ' num2str(num_in_b(3),'%2.1e') ],'Location','northwest')
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);



figure(ifig)
ifig=ifig+1;
cc=jet(length(aaa))
for ni=1:1:qq
hold on
plot(eval(['x_coor' int2str(ni) ]),eval(['y_coor' int2str(ni) ]),...
  '.','color',cc(ni,:),'LineWidth',0.5)

% plot(linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))),...
%     histc(eval(['full_spectrum_' int2str(ni) ]),linspace(0,max(eval(['full_spectrum_' int2str(ni) ]))))*weigth,...
%     'LineWidth',2);

grid on
hold off
end
set(gca,'FontSize',16)
ylabel('x')
xlabel('y')
legend('250nC','275nC','225nC',0)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);

































