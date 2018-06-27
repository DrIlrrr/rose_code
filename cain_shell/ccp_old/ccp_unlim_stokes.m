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


   STOKES(1,:)=[1 0 0];
   STOKES(2,:)=[0 1 0];
   STOKES(3,:)=[0 0 1];
    
    aaa=[1 2 3]
for var_for_scan=aaa;
    qq=qq+1;
 
    
    rflags.STOKES=STOKES(var_for_scan,:);
    

    
    
    rflags.pulseE=0.4; %laser puse energy [J]
    
    
    
    defoc_param=0.9;%defocusing parameter by defoult 1 is no defocusing
% BASE_DIRECTORY = [pwd '/dif_STOKES_' num2str(rflags.STOKES(1)) '_' num2str(rflags.STOKES(2)) '_' num2str(rflags.STOKES(3)) '_polar_pulseE_' num2str(rflags.pulseE) '_hi_600_defocus_' num2str(defoc_param) '_PulseE_400mJ/'];
        BASE_DIRECTORY = [pwd '/dif_STOKES_' num2str(rflags.STOKES(1)) '_' num2str(rflags.STOKES(2)) '_' num2str(rflags.STOKES(3))...
        '_polar_pulseE_' num2str(rflags.pulseE) '_hi_600_MeV_defocus_' num2str(defoc_param) '_PulseE_400mJ/'];
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
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);






theta=7.8e-5;
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
 title({['Spectrum of scattered photons in Theta=' num2str(theta) ' [rad]']})%,...
%     ['bandwith =' num2str(std(full_spectrum_1(aa_1))/mean(full_spectrum_1(aa_1))) ' for \delta t=0']})
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')
% legend(['\alpha_0=0 ; photons = ' num2str(num_in_b(1)*weigth,'%2.1e') ],...
%     ['\alpha_0=2.5; photons = ' num2str(num_in_b(2)*weigth,'%2.1e') ],...
%     ['\alpha_0=5; photons = ' num2str(num_in_b(3)*weigth,'%2.1e') ],...
%     ['\alpha_0=7.5; photons = ' num2str(num_in_b(4)*weigth,'%2.1e') ],0)
filename = ['plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);










figure(ifig)
ifig=ifig+1;
cc=jet(length(aaa))
for ni=1:1:3
    
    subplot(1,3,ni)
    
hold on
plot(eval(['x_coor' int2str(ni) ]),eval(['y_coor' int2str(ni) ]),...
  '.','LineWidth',0.1)

grid on
hold off
% xlim([-5e-5 5e-5])
% ylim([-5e-5 5e-5])
set(gca,'FontSize',16)
ylabel('x')
xlabel('y')

end




figure(ifig)
ifig=ifig+1;
for ni=1:1:3
xedges=[];
xedges = linspace(-1.5e-4,1.5e-4,1e2); yedges =linspace(-1.5e-4,1.5e-4,1e2);
histmat = hist2(eval(['x_coor' int2str(ni) ]),eval(['y_coor' int2str(ni) ]), xedges, yedges);


subplot(1,3,ni)
set(gca,'FontSize',16)
hold on
% mesh(xedges,yedges,histmat')
 set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
xlabel('x')
ylabel('y')
end




figure(ifig)
ifig=ifig+1;
for ni=1:1:3
xedges=[];
xedges = linspace(-0.6e4,0.6e4,1e2); yedges =linspace(-0.6e4,0.6e4,1e2);
histmat = hist2(eval(['xp_' int2str(ni) ]),eval(['yp_' int2str(ni) ]), xedges, yedges);

subplot(1,3,ni)
set(gca,'FontSize',16)
hold on
%  mesh(xedges,yedges,histmat')
 set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
xlabel('Px')
ylabel('Py')
end







l=1;

for ni=1:1:3

eval(['x1_coor' int2str(ni) '=x_coor' int2str(ni) '+(xp_' int2str(ni) './zp_' int2str(ni) ').*l;']);
eval(['y1_coor' int2str(ni) '=y_coor' int2str(ni) '+(yp_' int2str(ni) './zp_' int2str(ni) ').*l;']);
end

figure(ifig)
ifig=ifig+1;

for ni=1:1:3
xedges=[];
xedges = linspace(-40e-4,40e-4,1e2); yedges =linspace(-40e-4,40e-4,1e2);
histmat = hist2(eval(['x1_coor' int2str(ni) ]),eval(['y1_coor' int2str(ni) ]), xedges, yedges);

subplot(1,3,ni)
set(gca,'FontSize',16)
hold on
% mesh(xedges,yedges,histmat')
set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
xlabel('x')
ylabel('y')
suptitle(['after ' num2str(l) ' m'])
end



figure(ifig)
ifig=ifig+1;


for ni=1:1:3
xedges=[];
xedges = linspace(-0.4e-4,0.4e-4,1e2); yedges =linspace(-1e4,1e4,1e2);
histmat = hist2(eval(['x_coor' int2str(ni) ]),eval(['xp_' int2str(ni) ]), xedges, yedges);

subplot(1,3,ni)
set(gca,'FontSize',16)
hold on
% mesh(xedges,yedges,histmat')
set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
xlabel('x')
ylabel('Px')
% suptitle(['after ' num2str(l) ' m'])
end


figure(ifig)
ifig=ifig+1;


for ni=1:1:3
xedges=[];
xedges = linspace(-0.4e-4,0.4e-4,1e2); yedges =linspace(-1e4,1e4,1e2);
histmat = hist2(eval(['y_coor' int2str(ni) ]),eval(['yp_' int2str(ni) ]), xedges, yedges);

subplot(1,3,ni)
set(gca,'FontSize',16)
hold on
% mesh(xedges,yedges,histmat')
set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
xlabel('y')
ylabel('Py')
% suptitle(['after ' num2str(l) ' m'])
end





figure(ifig)
ifig=ifig+1;

for ni=1:1:3
    
xedges=[];
xedges = linspace(-1e-4,1e-4,1e2); yedges =linspace(-6e3,6e3,1e2);
histmat = hist2(eval(['x1_coor' int2str(ni) ]),eval(['xp_' int2str(ni) ]), xedges, yedges);

subplot(1,3,ni)
set(gca,'FontSize',16)
hold on
% mesh(xedges,yedges,histmat')
set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
xlabel('x')
ylabel('Px')
suptitle(['after ' num2str(l) ' m'])
end
set(gca,'FontSize',16)
ylabel('x')
xlabel('y')











%% Sx S y Sz


figure(ifig)
ifig=ifig+1;
for ni=1:1:3


subplot(1,3,ni)
hist(eval(['Sx_' int2str(ni) ]), 100);
xlabel('Sx')
% suptitle(['after ' num2str(l) ' m'])
end



figure(ifig)
ifig=ifig+1;
for ni=1:1:3


subplot(1,3,ni)
hist(eval(['Sy_' int2str(ni) ]), 100);
xlabel('Sy')
% suptitle(['after ' num2str(l) ' m'])
end



figure(ifig)
ifig=ifig+1;
for ni=1:1:3


subplot(1,3,ni)
hist(eval(['Sz_' int2str(ni) ]), 100);
xlabel('Sz')
% suptitle(['after ' num2str(l) ' m'])
end




figure(ifig)
ifig=ifig+1;

for ni=1:1:3
    
xedges=[];
xedges = linspace(-1,1,1e2); yedges =linspace(-1,1,1e2);
histmat = hist2(eval(['Sx_' int2str(ni) ]),eval(['Sz_' int2str(ni) ]), xedges, yedges);

subplot(1,3,ni)
set(gca,'FontSize',16)
hold on
% mesh(xedges,yedges,histmat')
set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
set(gca,'FontSize',16)
ylabel('Sx')
xlabel('Sz')
end



figure(ifig)
ifig=ifig+1;

for ni=1:1:3
    
xedges=[];
xedges = linspace(-1,1,1e2); yedges =linspace(-1,1,1e2);
histmat = hist2(eval(['Sx_' int2str(ni) ]),eval(['Sy_' int2str(ni) ]), xedges, yedges);

subplot(1,3,ni)
set(gca,'FontSize',16)
hold on
% mesh(xedges,yedges,histmat')
set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
set(gca,'FontSize',16)
ylabel('Sx')
xlabel('Sy')
end



figure(ifig)
ifig=ifig+1;

for ni=1:1:3
    
xedges=[];
xedges = linspace(-1,1,1e2); yedges =linspace(-1,1,1e2);
histmat = hist2(eval(['Sy_' int2str(ni) ]),eval(['Sz_' int2str(ni) ]), xedges, yedges);

subplot(1,3,ni)
set(gca,'FontSize',16)
hold on
% mesh(xedges,yedges,histmat')
set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
set(gca,'FontSize',16)
ylabel('Sy')
xlabel('Sz')
end



kuku=1:1:length(Sx_1);
figure(ifig)
ifig=ifig+1;
hold on
plot3(Sx_1(kuku),Sy_1(kuku),Sz_1(kuku),'.b')
 plot3(Sx_2(kuku),Sy_2(kuku),Sz_2(kuku),'.r')
plot3(Sx_3(kuku),Sy_3(kuku),Sz_3(kuku),'.c')
hold off
grid on
axis square



