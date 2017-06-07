function [full_spectrum,phot_angle,weigth,x_phot,y_phot,z_phot,xp_phot,yp_phot,zp_phot,Sx,Sy,Sz]=scan_photons_plots(out_folder)
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
zp_phot=photons_data(:,11);
Sx=photons_data(:,12);
Sy=photons_data(:,13);
Sz=photons_data(:,14);

weigth=max(photons_data(:,3))

full_spectrum=photons_data(:,8)./1e3;
% phot_angle=[phot_angle;sign(photons_data(:,9)).*atan(sqrt(photons_data(:,9).^2+photons_data(:,10).^2)./photons_data(:,11))];
phot_angle=sign(photons_data(:,9)).*atan(sqrt(photons_data(:,9).^2+photons_data(:,10).^2)./photons_data(:,11));

number_of_photons=length(find(full_spectrum))*weigth;
%SPEED_OF_LIGHT=3e8;


end