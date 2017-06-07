clear all; close all; clc;

%% Constants
SPEED_OF_LIGHT=3e8;
Sigma_th=0.665e-28;%[m^2] Thomsom cross section
h=2*pi*1.054e-34;%Planc const [J*s]
Echarge = 1.60e-19;% Charge of electron [c]
%% Laser parameteres
angle=0;
sigLr=30e-6;
sigt=6e-12;
lambda_l=800e-9;% [m] laser wave length
pulseE=5; %[J]
photons_number=pulseE/((h*SPEED_OF_LIGHT)/lambda_l);
laser_length=sigt*SPEED_OF_LIGHT;
%% electron bunch parameters
chargebunch = 250e-12;%Charge per electrons bunch [c]
std_x_size=[14*1e-6]; % IP vertical electron beam size [m];
std_y_size=std_x_size; % IP horizontal electron beam size [m]
std_s=880*1e-6;%6*1e-12*SPEED_OF_LIGHT; % intial bunch length [s]
number_electrons =chargebunch/Echarge;% Number electrons in bunch
%%
delta_s=0;
% delta_x=0;
delta_y=0;

figure(1)
hold on
qq=0;
delta_xaar=[0];
for qq=[1:1:length(delta_xaar)]
    delta_x=delta_xaar(qq);
Ax(qq,:)=exp((delta_x.^2)./(-2.*(std_x_size.^2+sigLr.^2+((std_s).^2+(sigt*SPEED_OF_LIGHT)^2).*tan((1/2)*angle*(pi/180))^2)));

Ay(qq,:)=exp((delta_y.^2)./(-2*(std_x_size.^2+sigLr^2)));

As(qq,:)=exp((delta_s.^2.*tan((1/2)*angle.*(pi/180))^2)./(-2*(std_x_size.^2+sigLr^2+((std_s)^2+(sigt*SPEED_OF_LIGHT)^2)*tan((1/2)*angle*(pi/180))^2)));

NUM_ph_x(qq,:)=Ax(qq,:).*Ay(qq,:).*As(qq,:).*((Sigma_th.*photons_number.*number_electrons)./(2.*pi.*sqrt(std_y_size.^2+sigLr.^2)))...
    ./sqrt((std_x_size.^2+sigLr.^2)+((std_s).^2+(sigt.*SPEED_OF_LIGHT).^2).*tan((1/2)*angle*(pi/180))^2);

plot(std_x_size*1e6,NUM_ph_x(qq,:),'*r')

end
hold off
grid on 
xlabel('\sigma_x [\mu m]','FontSize',16)
ylabel('Flux','FontSize',16)
 legendCell = cellstr(num2str(delta_xaar', 'dx=%-10.2e'))
legend(legendCell,'FontSize',16)
