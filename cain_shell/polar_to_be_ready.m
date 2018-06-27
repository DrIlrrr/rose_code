%% for analisys of polarisation

clear all; close all; clc;

%make_path
start_date=datestr(now);
%for new global
global home_dir DIRECTORY_FOR_CAIN BASE_DIRECTORY;
global rflags
[rflags] = flags_for_run;
out_folder=[pwd '/polar_plot/'];
mkdir(out_folder)
ifig=1;
rflags.PLOTS =1;
just_plots=1;
qq=0;
scan_spectrum=[];
phot_angle=[];


%STOKES(1,:)=[1 0 0];
%STOKES(2,:)=[0 1 0];
%STOKES(3,:)=[0 0 1];

rflags.pulseE=0.2; %laser puse energy [J]

rflags.chargebunch = 250e-12;%Charge per electrons bunch [c]
defoc_param=1;%defocusing parameter by defoult 1 is no defocusing
BASE_DIRECTORY = [pwd '/'];% '/Echarge_' num2str(rflags.chargebunch*1e12) '_pulseE_' num2str(rflags.pulseE) '/'];
%beam_phasespace=dlmread(['cristina_e_beams/eli_highen_oned_WP_newlayout_track_up_new_check_newsol_rec_600MeV.out.asci'],'',10,0);

%%%
[full_spectrum,phot_angle,weigth,x_phot,y_phot,z_phot,xp_phot,yp_phot,zp_phot,Sx,Sy,Sz]=scan_photons_plots(BASE_DIRECTORY);

%% after propagation


bu_factor=100
theta=0.000192;
w_size=bu_factor*5e-5;
nx_bin=50;
ny_bin=50;
xcor_allb=[]; ycor_allb=[]; x1_coor=[]; y1_coor=[];
aa=find(abs(phot_angle)<theta);


for l=1;%[0:0.1:1]
    close all
    
    xcor_allb=x_phot+(xp_phot./zp_phot).*l;
    ycor_allb=y_phot+(yp_phot./zp_phot).*l;
    
    Sz_sep=Sz(aa);
    
    xcor_n=x_phot(aa)+(xp_phot(aa)./zp_phot(aa)).*l;
    ycor_n=y_phot(aa)+(yp_phot(aa)./zp_phot(aa)).*l;
    
    d1=zeros(ny_bin,nx_bin);
    d2=zeros(ny_bin,nx_bin);
    f=zeros(ny_bin,nx_bin);
    d1_theta=zeros(ny_bin,nx_bin);
    d2_theta=zeros(ny_bin,nx_bin);
    f_theta=zeros(ny_bin,nx_bin);
    
    x_l=2*w_size/nx_bin;
    y_l=2*w_size/ny_bin;
    
    ang=0:0.01:2*pi;
    r=l*sin(theta)
    xpt=r*cos(ang);
    ypt=r*sin(ang);
    
    %%
    xedges = linspace(-5e-4,5e-4,1e2); yedges =linspace(-5e-4,5e-4,1e2);
    histmat = hist2(xcor_allb,ycor_allb, xedges, yedges);
    
    figure(ifig)
    ifig=ifig+1;
    set(gca,'FontSize',16)
    % mesh(xedges,yedges,histmat')
    set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    hold on
    plot(xpt,ypt,'-k')
    hold off
    colormap(jet)
    xlabel('x')
    ylabel('y')
    title(['at L=' num2str(l) ' m'])
    filename = [ out_folder  'XY_plot_' num2str(ifig) ];
    fname = [ filename '.png'];
    print('-dpng', fname);
    
    
    
    %%
    
    
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
    
    
    ang=0:0.01:2*pi;
    r=l*sin(theta)
    xpt=r*cos(ang);
    ypt=r*sin(ang);
    
    figure(100)
    % ifig=ifig+1;
    subplot 121
    set(pcolor(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        f),'EdgeColor','none')
    colormap(jet)
    hold on
    plot(xpt,ypt,'-r')
    hold off
    colorbar('SouthOutside')
    title('prop')
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    subplot 122
    set(pcolor(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        d2),'EdgeColor','none')
    hold on
    plot(xpt,ypt,'-r')
    hold off
    colorbar('SouthOutside')
    title('prop')
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    suptitle(['L=' num2str(l) ' m; Stokes=' num2str(rflags.STOKES) ])
    filename = [ out_folder 'EXP_plot_' num2str(1) '_' num2str(l*1e6) ];
    fname = [ filename '.png'];
    print('-dpng', fname);
    
    
    
    
    figure(200)
    % ifig=ifig+1;
    subplot 221%121
    set(pcolor(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        f),'EdgeColor','none')
    colormap(jet)
    hold on
    plot(xpt,ypt,'-r')
    hold off
    colorbar('SouthOutside')
    title('prop')
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    
    subplot 222%122
    set(pcolor(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        d2),'EdgeColor','none')
    hold on
    plot(xpt,ypt,'-r')
    hold off
    colorbar('SouthOutside')
    title('prop')
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    %     suptitle(['L=' num2str(l) ' m; Stokes=' num2str(rflags.STOKES) ])
    
    %         figure(200)
    subplot 223%121
    
    plot(0,0,'xr')
    set(pcolor(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        f_theta),'EdgeColor','none')
    hold on
    plot(xpt,ypt,'-r')
    hold off
    colorbar('SouthOutside')
    title('prop')
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    
    subplot 224%122
    set(pcolor(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        d2_theta),'EdgeColor','none')
    hold on
    plot(xpt,ypt,'-r')
    hold off
    colorbar('SouthOutside')
    title('prop')
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    suptitle(['L=' num2str(l) ' m; Stokes=' num2str(rflags.STOKES) ])
    
    
    filename = [ out_folder 'EXP_plot_' num2str(2) '_' num2str(l*1e6) ];
    fname = [ filename '.png'];
    print('-dpng', fname);
    
end







xedges = linspace(-1.5e-4,1.5e-4,1e2); yedges =linspace(-1.5e-4,1.5e-4,1e2);
histmat = hist2(x_phot,y_phot, xedges, yedges);

figure(ifig)
ifig=ifig+1;
set(gca,'FontSize',16)
hold on
% mesh(xedges,yedges,histmat')
set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
colormap(jet)
xlabel('x')
ylabel('y')
title('at IP')
filename = [ out_folder  'XY_IP_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);








bandwith_cm=[];
num_phot_th=[];
theta_angle=10e-6;
diapason=(10:1:101);
el_angel=0;
qq=0; num_phot_th=[]; bandwith_cm=[]; bandwith_non_norm=[];
for ni=diapason;
    qq=qq+1;
    el_angel_1=ni*theta_angle;
    aa_1=find(abs(phot_angle)<el_angel_1);
    
    num_phot_th(qq)=length(full_spectrum(aa_1))*weigth;
    bandwith_cm(qq)=std(full_spectrum(aa_1))/mean(full_spectrum(aa_1));
    
    mean_Sx(qq)=mean(Sx(aa_1));
    mean_Sy(qq)=mean(Sy(aa_1));
    mean_Sz(qq)=mean(Sz(aa_1));
    
    if bandwith_cm(qq)>=0.0048 && bandwith_cm(qq)<=0.0051 % find angle for given bandwith
        el_angel=el_angel_1
        % stop;
    end
end
aa=find(abs(phot_angle)<el_angel);

bandwith=std(full_spectrum(aa))/mean(full_spectrum(aa))
el_angel



figure(ifig)
ifig=ifig+1;
set(gca,'FontSize',16)
subplot(3,3,[1 3])
plot(diapason.*theta_angle,num_phot_th,'-xb')
grid on
set(gca,'FontSize',16)
ylabel('number photons ')
subplot(3,3,[4 6])
hold on
plot(diapason.*theta_angle,bandwith_cm,'--xb')
plot(el_angel,bandwith,'or','LineWidth',2)
hold off
set(gca,'FontSize',16)
ylabel('bandwidth')
xlabel('Theta')
grid on
set(gca,'FontSize',16)
subplot(3,3,7)
plot(diapason.*theta_angle,mean_Sx,'-xb')
grid on
ylabel('mean Sx')
xlabel('Theta')
 subplot(3,3,8)
 plot(diapason.*theta_angle,mean_Sy,'-xb')
 grid on
ylabel('mean Sy')
xlabel('Theta')
 subplot(3,3,9)
 plot(diapason.*theta_angle,mean_Sz,'-xb')
 grid on
ylabel('mean Sx')
xlabel('Theta')

suptitle({['bandwith=' num2str(bandwith) ' for \theta=' num2str(el_angel) ];...
    ['full flux=' num2str(length(find(full_spectrum))) ';'...
    ' Flux in \theta =' num2str(length(find(full_spectrum(aa))),'%10.2e') ]})
filename = [ out_folder  'photons_plot_' num2str(ifig,'%10.2e') ];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
print('-depsc', fname);print('-dpng', fname2);



figure(ifig)
ifig=ifig+1;
set(gca,'FontSize',16)
subplot(5,1,1)
plot(diapason.*theta_angle,num_phot_th,'-xb')
grid on
set(gca,'FontSize',16)
ylabel('flux')
subplot(5,1,2)
hold on
plot(diapason.*theta_angle,bandwith_cm,'--xb')
plot(el_angel,bandwith,'or','LineWidth',2)
hold off
set(gca,'FontSize',16)
ylabel('BD')
xlabel('Theta')
grid on
set(gca,'FontSize',16)
subplot(5,1,3)
plot(diapason.*theta_angle,mean_Sx,'-xb')
grid on
ylabel('mean Sx')
xlabel('Theta')
set(gca,'FontSize',16)
subplot(5,1,4)
 plot(diapason.*theta_angle,mean_Sy,'-xb')
 grid on
ylabel('mean Sy')
xlabel('Theta')
set(gca,'FontSize',16)
subplot(5,1,5)
 plot(diapason.*theta_angle,mean_Sz,'-xb')
 grid on
ylabel('mean Sx')
xlabel('Theta')
set(gca,'FontSize',16)
suptitle({['bandwith=' num2str(bandwith) ' for \theta=' num2str(el_angel) ];...
    ['full flux=' num2str(length(find(full_spectrum))) ';'...
    ' Flux in \theta =' num2str(length(find(full_spectrum(aa))),'%10.2e') ]})
filename = [ out_folder  'photons_plot_' num2str(ifig,'%10.2e') ];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
print('-depsc', fname);print('-r300','-dpng', fname2);

 figure(1200)
    set(pcolor(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        f),'EdgeColor','none')
    colormap(jet)
    hold on
    plot(xpt,ypt,'-r','LineWidth',2)
    hold off
    colorbar('SouthOutside')
     title(['L=' num2str(l) ' m; Stokes=' num2str(rflags.STOKES) ';'])
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    filename = [ out_folder  'S3_col' ];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
print('-depsc', fname);print('-dpng', fname2);
    
    
    
  figure(1201)
    set(pcolor(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        d2),'EdgeColor','none')
    colormap(jet)
    hold on
    plot(xpt,ypt,'-k','LineWidth',2)
    hold off
    colorbar('SouthOutside')
title(['L=' num2str(l) ' m; Stokes=' num2str(rflags.STOKES) ';'])
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    %suptitle(['L=' num2str(l) ' m; Stokes=' num2str(rflags.STOKES) ])
        filename = [ out_folder  'XY_col' ];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
print('-depsc', fname);print('-dpng', fname2);





figure(1300)
mesh(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        f)
     colormap(jet)
    hold on
    plot(xpt,ypt,'-r','LineWidth',2)
    hold off
    colorbar('SouthOutside')
     title(['L=' num2str(l) ' m; Stokes=' num2str(rflags.STOKES) ';'])
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    filename = [ out_folder  'S3_col' ];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
print('-depsc', fname);print('-dpng', fname2);





figure(ifig)
ifig=ifig+1;
set(gca,'FontSize',16)
subplot(3,1,1)
plot(diapason.*theta_angle,num_phot_th,'-xb')
grid on
set(gca,'FontSize',16)
ylabel('flux')
subplot(3,1,2)
hold on
plot(diapason.*theta_angle,bandwith_cm,'--xb')
plot(el_angel,bandwith,'or','LineWidth',2)
hold off
set(gca,'FontSize',16)
ylabel('BD')
% xlabel('Theta')
grid on
set(gca,'FontSize',16)
set(gca,'FontSize',16)
subplot(3,1,3)
 plot(diapason.*theta_angle,mean_Sz,'-xb')
 grid on
ylabel('mean Sz')
xlabel('Theta')
set(gca,'FontSize',16)
suptitle({['bandwith=' num2str(bandwith) ' for \theta=' num2str(el_angel) ];...
    ['full flux=' num2str(length(find(full_spectrum))) ';'...
    ' Flux in \theta =' num2str(length(find(full_spectrum(aa))),'%10.2e') ]})
filename = [ out_folder  'photons_plot_' num2str(ifig,'%10.2e') ];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
print('-depsc', fname);print('-r300','-dpng', fname2);

%%%
nbin_plot=500;
xs=linspace(min(xcor_n),max(xcor_n),nbin_plot);
ys=smooth(hist(xcor_n,nbin_plot));
iysm=find(ys==max(ys));

figure(ifig)
ifig=ifig+1;
set(gca,'FontSize',16)
hold on
hist(xcor_n,nbin_plot)
plot(xs,ys,'-r','LineWidth',3)
plot(xs(iysm/2),ys(iysm/2),'or','LineWidth',5)
hold off
grid on
xlabel('X [m]')
title({['L=' num2str(l) ' m; Stokes=' num2str(rflags.STOKES) ';'];[' BW=' num2str(bandwith) ' for \theta=' num2str(el_angel) ]})

