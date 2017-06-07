clear all; close all; clc;
% electron_data=load([pwd '/try_long_laser_shifting_0/cain_tmp/exp.dat']);
std_x_size=17.2*1e-6;%std(electron_data(:,5));
std_y_size=16.4*1e-6;%std(electron_data(:,6));
std_s=2.7e-4;%std(electron_data(:,7));
Echarge = 1.60e-19;% Charge of electron [c]
number_electrons = 250e-12/Echarge;% Number electrons in bunch

SPEED_OF_LIGHT=3e8;
Sigma_th=0.665e-28;%[m^2] Thomsom cross section
h=2*pi*1.054e-34;%Planc const [J*s]


lambda_l=515e-9;% [m] laser wave length
pulseE=0.4; %[J]
photons_number=pulseE/((h*SPEED_OF_LIGHT)/lambda_l);

angle=0;
sigLr=28e-6/2;
sigt=1.5e-12;
% delta_s=0;
% delta_x=0;
% 
% As=exp((delta_s^2.*tan((1/2)*angle.*(pi/180)))/(-2*((std_x_size.^2+sigLr^2)+((std_s).^2+(sigt.*SPEED_OF_LIGHT).^2).*tan((1/2)*angle.*(pi/180)))));
% Ax=exp((delta_x^2)/(-2*((std_x_size.^2+sigLr^2)+((std_s).^2+(sigt.*SPEED_OF_LIGHT).^2).*tan((1/2)*angle.*(pi/180)))));



NUM_ph=((Sigma_th*photons_number*number_electrons)/(2*pi*sqrt(std_y_size.^2+sigLr^2)))...
    /sqrt((std_x_size.^2+sigLr^2)+((std_s).^2+(sigt.*SPEED_OF_LIGHT).^2).*tan((1/2)*angle.*(pi/180)))




