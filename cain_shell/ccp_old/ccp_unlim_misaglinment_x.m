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


aaa=[-100:20:100];


for var_for_scan=aaa
    qq=qq+1;
    
    rflags.shifting_laser_x = var_for_scan*1e-6;  %shifting laser possition at IP [m]
    
    rflags.pulseE=0.2; %laser puse energy [J]
    
    rflags.chargebunch = 250e-12;%Charge per electrons bunch [c]
    defoc_param=1;%defocusing parameter by defoult 1 is no defocusing
    
    BASE_DIRECTORY = [pwd '/Echarge_' num2str(rflags.chargebunch*1e12) '_pulseE_' num2str(rflags.pulseE) '_x_shift_' num2str(rflags.shifting_laser_x) '/'];
    
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
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    electron_data=load([BASE_DIRECTORY 'cain_tmp/exp.dat']);
    std_x_size=std(electron_data(:,5));
    std_y_size=std(electron_data(:,6));
    std_s=std(electron_data(:,4));
    number_electrons=electron_data(1,3)*length(electron_data(:,3));
    
    
    std_x_size=std(electron_data(:,5));
    std_y_size=std(electron_data(:,6));
    std_s=std(electron_data(:,4));
    number_electrons=electron_data(1,3)*length(electron_data(:,3));
    
    SPEED_OF_LIGHT=3e8;
    Sigma_th=0.665e-28;%[m^2] Thomsom cross section
    h=2*pi*1.054e-34;%Planc const [J*s]
    lambda_l=515e-9;% [m] laser wave length
    pulseE=0.2; %[J]
    photons_number=pulseE/((h*SPEED_OF_LIGHT)/lambda_l);
    
    angle=8;
    sigLr=14e-6;
    %         laser_length=var_val;
    sigt=1.5e-12;
    laser_length=sigt*SPEED_OF_LIGHT;
    
    delta_s=0;
    delta_x(qq)=rflags.shifting_laser_x;
%     Ax(qq)=exp((delta_x(qq).^2)/(-2*(std_y_size.^2+sigLr^2)));
    
     Ax(qq)=exp((delta_x(qq).^2)/(-2*(std_x_size.^2+sigLr^2+((std_s)^2+(sigt*SPEED_OF_LIGHT)^2)*tan((1/2)*angle*(pi/180))^2)));
    NUM_ph_x(qq)=Ax(qq).*((Sigma_th*photons_number*number_electrons)/(2*pi*sqrt(std_y_size^2+sigLr^2)))...
        /sqrt((std_x_size^2+sigLr^2)+((std_s)^2+(sigt*SPEED_OF_LIGHT)^2)*tan((1/2)*angle*(pi/180))^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
end


color_line={'-r','-b','-g','-m','-y','--y','-r','-b','-g','-m','-y','--y','-r','-b','-g','-m','-y','--y','-r','-b','-g','-m','-y','--y'};
color_dot={'-r','-b','-g','-m','-y','--y','-r','-b','-g','-m','-y','--y','-r','-b','-g','-m','-y','--y','-r','-b','-g','-m','-y','--y'};

figure(ifig)
ifig=ifig+1;
plot(delta_x,NUM_ph_x,'LineWidth',2)
xlabel('\Delta x','FontSize',16)
ylabel('Flux','FontSize',16)
grid on
% filename = [ 'errors_' num2str(ifig) ];


figure(ifig)
ifig=ifig+1;
for ni=1:1:qq
    qq
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
legend('0 \mu m','20 \mu m','40 \mu m','60 \mu m','80 \mu m','100 \mu m',0)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);


ww=0;
for ni=aaa
    ww=ww+1
    eval(['aa_' num2str(ww) '=find(abs(phot_angle_' num2str(ww) '));']);
    total_num_in_b(ww)=eval(['length(full_spectrum_' int2str(ww) '(aa_' int2str(ww) '))*weigth_' int2str(ww) ';']);
    
end

theta=19.5e-5;
ww=0;
for ni=aaa
    ww=ww+1
    eval(['aa_' num2str(ww) '=find(abs(phot_angle_' num2str(ww) ')<' num2str(theta) ');']);
    num_in_b(ww)=eval(['length(full_spectrum_' int2str(ww) '(aa_' int2str(ww) '))*weigth_' int2str(ww) ';']);
    
end



figure(ifig)
ifig=ifig+1;
hold on
plot(delta_x,NUM_ph_x,'LineWidth',2)
plot(delta_x,total_num_in_b,'-.xr','LineWidth',0.5)
plot(delta_x,total_num_in_b,'xr','LineWidth',3)
hold off
xlabel('\Delta x','FontSize',16)
ylabel('Flux','FontSize',16)
grid on
set(gca,'FontSize',16)
legend('theor','cain','Location','NorthEast')
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);

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
ylabel('Flux')%('number of scattered photons')
subplot 212
hold on
plot(aaa,num_in_b,'-.xb','LineWidth',0.5)
plot(aaa,num_in_b,'xb','LineWidth',3)
title(['in Theta=' num2str(theta) ])
hold off
grid on
ylim([0 max(num_in_b)+max(num_in_b)*1e-1])
set(gca,'FontSize',16)
ylabel('Flux')%('number of scattered photons')
xlabel('X shifting \mu m')
% legend('\delta t=0 [ps]','\delta t=1 [ps]','\delta t=2 [ps]',0)
filename = ['plot_' num2str(ifig) ];% one beam several events
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
legend(['0; photons = ' num2str(num_in_b(1),'%2.1e') ],...
    ['20 \mu m; photons = ' num2str(num_in_b(2),'%2.1e') ],...
    ['40 \mu m; photons = ' num2str(num_in_b(3),'%2.1e') ],...
    ['60 \mu m; photons = ' num2str(num_in_b(4),'%2.1e') ],...
    ['80 \mu m; photons = ' num2str(num_in_b(5),'%2.1e') ],...
    ['100 \mu m; photons = ' num2str(num_in_b(6),'%2.1e') ],'Location','northwest')
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);



vec1=[delta_x; total_num_in_b; num_in_b];
f_name=['mis_x.txt'];
fileID = fopen(f_name,'w');
fprintf(fileID,'%10.5e %10.5e %10.5e\n',vec1);
fclose(fileID);

