clear all; close all; clc;
make_path

home_dir=[pwd '/CAIN/'];
mkdir(home_dir)

indic=0;
turn_number=20
for aa=[0 1 2 3 4 5 6 7 7.5]
indic=indic+1;

% out_folder = [pwd '/try_long_laser_shifting_' num2str(rflags.shifting_laser_t) '/'];
out_folder = [pwd '/new_beam_' num2str(aa) '/'];

electron_data=load([out_folder 'cain_tmp/exp.dat']);
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

angle=aa;
sigLr=28e-6/2.0;
%         laser_length=var_val;
sigt=1.5e-12;
laser_length=sigt*SPEED_OF_LIGHT;


% As=exp((delta_s^2.*tan((1/2)*angle.*(pi/180)))/(-2*((std_x_size.^2+sigLr^2)+((std_s).^2+(sigt.*SPEED_OF_LIGHT).^2).*tan((1/2)*angle.*(pi/180)))));
% Ax=exp((delta_x^2)/(-2*((std_x_size.^2+sigLr^2)+((std_s).^2+(sigt.*SPEED_OF_LIGHT).^2).*tan((1/2)*angle.*(pi/180)))));
%
%

% 
NUM_ph=((Sigma_th*photons_number*number_electrons)/(2*pi*sqrt(std_y_size^2+sigLr^2)))...
    /sqrt((std_x_size^2+sigLr^2)+((std_s)^2+(sigt*SPEED_OF_LIGHT)^2)*tan((1/2)*angle*(pi/180))^2);
% 
% NUM_ph=((Sigma_th*photons_number*number_electrons)/(2*pi*sqrt(std_x_size^2+sigLr^2)))...
%     /sqrt((std_y_size^2+sigLr^2)+((std_s)^2+(sigt*SPEED_OF_LIGHT)^2)*tan((angle/2)));


all_data_phot=[];
full_spectrum=[];
phot_angle=[];
weigth=0;
histmat=[];
new_m=[];
x_phot=[];
y_phot=[];
xp_phot=[];
yp_phot=[];

for ni=1:1:turn_number
    
    n_turn(ni)=NUM_ph;
    
    
    check_photons_output=dir([out_folder 'cain_tmp/cain_output_photons_' num2str(ni) '.dat']);% checking cain output for photons
    
    check_photons_output_size=check_photons_output(:).bytes;
    
    if (check_photons_output_size==0)%if cain non produce photons write photons_property=[0,0,0,0,0,0,0,0,0,0,0,0,0,0];
        photons_data=zeros(1,14);
    else
        
        photons_data=dlmread([out_folder 'cain_tmp/cain_output_photons_' num2str(ni) '.dat'],'',1,0);%read electrons from cain
    end
    
    
    
    weigth=photons_data(1,3);
    number_of_photons(ni)=length(photons_data)*weigth;
    
    
end

figure(indic)
hold on
plot(number_of_photons,'--xb')
plot(n_turn,'-or')
hold off
number_of_photons

NUM_ph

n_cain(indic)=mean(number_of_photons)
n_formula(indic)=NUM_ph
end


figure(10000)
hold on
plot([0 1 2 3 4 5 6 7 7.5],n_cain,'--xb')
plot([0 1 2 3 4 5 6 7 7.5],n_formula,'-or')
hold off
grid on
set(gca,'FontSize',16)
legend('cain','formula',0)
xlabel('interaction anlge')
ylabel('number of scattered photons')
fname = ['compare_flux.png'];
print('-dpng', fname);
