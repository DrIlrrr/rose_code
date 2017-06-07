
function [X_RAY_photons_property] = read_photons_data
%read photons data from cain output

global   DIRECTORY_FOR_CAIN;

check_exist=exist([ DIRECTORY_FOR_CAIN 'cain_output_photons.dat']);

if(check_exist==0)
X_RAY_photons_property=zeros(1,14);
else

check_photons_output=dir([DIRECTORY_FOR_CAIN 'cain_output_photons.dat']);% checking cain output for photons

check_photons_output_size=check_photons_output(:).bytes;

if (check_photons_output_size==0)%if cain non produce photons write photons_property=[0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    X_RAY_photons_property=zeros(1,14);
else
    
X_RAY_photons_property=dlmread([ DIRECTORY_FOR_CAIN 'cain_output_photons.dat'],'',1,0);


% system(['cat ' DIRECTORY_FOR_CAIN 'cain_output_photons.dat | sed -e 1d > ' DIRECTORY_FOR_CAIN 'cain_output_photons_' num2str(turn_number) '.dat']);% Delite first string (with text)
% X_RAY_photons_property=load([DIRECTORY_FOR_CAIN 'cain_output_photons_' num2str(turn_number) '.dat']); % Read file
end
end
%X_RAY_phot_energy=X_RAY_photons_property(:,ENERGY_OF_PARTICLE);
% 
% number_of_phot1=size(X_RAY_phot_energy);
% number_of_scatered_photons(turn_number)=number_of_phot1(1);
  
end

