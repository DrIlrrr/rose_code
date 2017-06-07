clear all; close all; clc;
electron_data=load([pwd '/try_long_laser_shifting_0/cain_tmp/exp.dat']);
std_x_size=std(electron_data(:,5));
std_y_size=std(electron_data(:,6));
std_s=std(electron_data(:,7));


SPEED_OF_LIGHT=3e8;
Sigma_th=0.665e-28;%[m^2] Thomsom cross section
h=2*pi*1.054e-34;%Planc const [J*s]
lambda_l=532e-9;% [m] laser wave length
pulseE=1; %[J]
photons_number=pulseE/((h*SPEED_OF_LIGHT)/lambda_l);

angle=0;
sigLr=30e-6/2.0;
%         laser_length=var_val;
laser_length=5e-12*SPEED_OF_LIGHT;
delta_s=rflags.shifting_laser_t;
delta_x=0;

As=exp((delta_s^2.*tan((1/2)*angle.*(pi/180)))/(-2*((std_x_size(index_num).^2+sigLr^2)+((std_s(index_num)).^2+(sigt.*SPEED_OF_LIGHT).^2).*tan((1/2)*angle.*(pi/180)))));
Ax=exp((delta_x^2)/(-2*((std_x_size(index_num).^2+sigLr^2)+((std_s(index_num)).^2+(sigt.*SPEED_OF_LIGHT).^2).*tan((1/2)*angle.*(pi/180)))));




NUM_ph(index_num)=As.*Ax.*((Sigma_th*photons_number*number_electrons)/(2*pi*sqrt(std_y_size(index_num).^2+sigLr^2)))...
    /sqrt((std_x_size(index_num).^2+sigLr^2)+((std_s(index_num)).^2+(sigt.*SPEED_OF_LIGHT).^2).*tan((1/2)*angle.*(pi/180)));




