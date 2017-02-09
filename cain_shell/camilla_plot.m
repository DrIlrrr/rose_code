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

aaa=[3]
for var_for_scan=aaa;
    qq=qq+1;
    
    
    rflags.STOKES=STOKES(var_for_scan,:);
    rflags.pulseE=0.4; %laser puse energy [J]
    
    beam_phasespace=dlmread(['cristina_e_beams/eli_highen_oned_WP_newlayout_track_up_new_check_newsol_rec_600MeV.out.asci'],'',10,0);
    beam_phasespace(:,6)=beam_phasespace(:,6)-114.7/0.511;
    defoc_param=0.9;%defocusing parameter by defoult 1 is no defocusing
    add_param=std(beam_phasespace(:,6));
    BASE_DIRECTORY = [pwd '/dif_STOKES_' num2str(rflags.STOKES(1)) '_' num2str(rflags.STOKES(2)) '_' num2str(rflags.STOKES(3))...
        '_polar_pulseE_' num2str(rflags.pulseE) '_hi_600_MeV_defocus_' num2str(defoc_param) '_PulseE_400mJ/'];
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


%% no propagation


w_size=5e-5;

nx_bin=100;
ny_bin=100;

d1=zeros(ny_bin,nx_bin);
d2=zeros(ny_bin,nx_bin);
f=zeros(ny_bin,nx_bin);


x_l=2*w_size/nx_bin;
y_l=2*w_size/ny_bin;


for ni=1:1:length(full_spectrum);
    if (x_coor1(ni)<w_size && x_coor1(ni)>-w_size)
        if (y_coor1(ni)<w_size && y_coor1(ni)>-w_size)
            d1(1+floor((y_coor1(ni)+w_size)./y_l),1+floor((x_coor1(ni)+w_size)./x_l))=d1(1+floor((y_coor1(ni)+w_size)./y_l),1+floor((x_coor1(ni)+w_size)./x_l))+Sz_1(ni);
            d2(1+floor((y_coor1(ni)+w_size)./y_l),1+floor((x_coor1(ni)+w_size)./x_l))=d2(1+floor((y_coor1(ni)+w_size)./y_l),1+floor((x_coor1(ni)+w_size)./x_l))+1;
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



% figure(ifig)
% ifig=ifig+1;
% %  pcolor(F)
% set(pcolor(d1),'EdgeColor','none')
% colorbar
%
%
% figure(ifig)
% ifig=ifig+1;
% %  pcolor(F)
% set(pcolor(d2),'EdgeColor','none')
% colorbar
%


figure(ifig)
ifig=ifig+1;
%  pcolor(F)
set(pcolor(linspace(-w_size,w_size,nx_bin),...
    linspace(-w_size,w_size,ny_bin),...
    f),'EdgeColor','none')
colorbar
title('IP')


figure(ifig)
ifig=ifig+1;
%  pcolor(F)
mesh(linspace(-w_size,w_size,nx_bin),...
    linspace(-w_size,w_size,ny_bin),...
    -f)
colorbar
title('IP')

%% after propagation



bu_factor=5
theta=bu_factor*11e-5;
ww=0;
for ni=3
    ww=ww+1
    eval(['aa_' num2str(ww) '=find(abs(phot_angle_' num2str(ww) ')<' num2str(theta) ');']);
    
end
ycor_nop=y_coor1(aa_1);
xcor_nop=x_coor1(aa_1);
Sz_sep=Sz_1(aa_1);

for l=1%[0:0.0001:0.01]% 0.1:0.1:0.5];%
    close all
x1_coor1=[];
y1_coor1=[];
for ni=1:1:1
    
    eval(['x1_coor' int2str(ni) '=x_coor' int2str(ni) '+(xp_' int2str(ni) './zp_' int2str(ni) ').*l;']);
    eval(['y1_coor' int2str(ni) '=y_coor' int2str(ni) '+(yp_' int2str(ni) './zp_' int2str(ni) ').*l;']);
    
        xcor=xcor_nop+(xp_1(aa_1)./zp_1(aa_1)).*l;
        ycor=ycor_nop+(yp_1(aa_1)./zp_1(aa_1)).*l;
    
    
end

w_size=bu_factor*11e-5;
nx_bin=50;
ny_bin=50;

d1=zeros(ny_bin,nx_bin);
d2=zeros(ny_bin,nx_bin);
f=zeros(ny_bin,nx_bin);
d1_theta=zeros(ny_bin,nx_bin);
d2_theta=zeros(ny_bin,nx_bin);
f_theta=zeros(ny_bin,nx_bin);

x_l=2*w_size/nx_bin;
y_l=2*w_size/ny_bin;


for ni=1:1:length(full_spectrum);
    if (x1_coor1(ni)<w_size && x1_coor1(ni)>-w_size)
        if (y1_coor1(ni)<w_size && y1_coor1(ni)>-w_size)
            d1(1+floor((y1_coor1(ni)+w_size)./y_l),1+floor((x1_coor1(ni)+w_size)./x_l))=...
                d1(1+floor((y1_coor1(ni)+w_size)./y_l),1+floor((x1_coor1(ni)+w_size)./x_l))+Sz_1(ni);
            d2(1+floor((y1_coor1(ni)+w_size)./y_l),1+floor((x1_coor1(ni)+w_size)./x_l))=...
                d2(1+floor((y1_coor1(ni)+w_size)./y_l),1+floor((x1_coor1(ni)+w_size)./x_l))+1;
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



for ni=1:1:length(xcor);
    if (xcor(ni)<w_size && xcor(ni)>-w_size)
        if (ycor(ni)<w_size && ycor(ni)>-w_size)
            d1_theta(1+floor((ycor(ni)+w_size)./y_l),1+floor((xcor(ni)+w_size)./x_l))=d1_theta(1+floor((ycor(ni)+w_size)./y_l),1+floor((xcor(ni)+w_size)./x_l))+Sz_sep(ni);
            d2_theta(1+floor((ycor(ni)+w_size)./y_l),1+floor((xcor(ni)+w_size)./x_l))=d2_theta(1+floor((ycor(ni)+w_size)./y_l),1+floor((xcor(ni)+w_size)./x_l))+1;
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
suptitle(['L=' num2str(l)])


filename = ['EXP_plot_' num2str(1) '_' num2str(l*1e6) ];
fname = [ filename '.png'];
print('-dpng', fname);

end





