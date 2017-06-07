function [beam_phasespace,X_RAY_photons_property] = start_cain_shell(beam_phasespace,turn_number)%,factor_def)


global  home_dir DIRECTORY_FOR_CAIN;
% global  beam_parameters
 global rflags

% plots_beam(beam_phasespace)

% beam_phasespace=dlmread(['eli_360MeV_double_A.out.asci'],'',10,0);
% a=factor_def;%defocusing parameter
% beam_phasespace(:,1)=beam_phasespace(:,1).*a;
% beam_phasespace(:,1)=beam_phasespace(:,1)./a;
% beam_phasespace(:,3)=beam_phasespace(:,3).*a;
% beam_phasespace(:,1)=beam_phasespace(:,1)./a;


% [beam_property]=formating_beam_for_cain(beam_phasespace);
[beam_property]=formating_beam_for_cain(beam_phasespace,turn_number);

[nothing] = start_cain(beam_property,turn_number);
X_RAY_photons_property=[];
%[X_RAY_photons_property] = read_photons_data;

%[beam_phasespace] = formating_beam_from_cain(beam_property);

% if(SAVE_ALL_DATA_FOR_PHOTONS==1)
%     save([FILE_BASENAME 'X_RAY_data_' num2str(turn_number) '_' num2str(NUMBER_OF_MACROPARTICLES) '_' num2str(pulseE) '.dat'],'X_RAY_photons_property','pulseE');
% end