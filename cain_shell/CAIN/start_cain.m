function [beam_property] = start_cain(beam_property,turn_number)
global  home_dir DIRECTORY_FOR_CAIN;
% global  beam_parameters
global rflags
DEBUGGING_MODE=0;

% number_of_steps_cain=250;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen([DIRECTORY_FOR_CAIN 'exp.dat'],'w');%save beam for cain standart
fprintf(fid,' %i    %i       %1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e \n',beam_property);
fclose(fid);
%%%%%%%%%%%%%%%% MODIFY INITIAL FILE FOR CAIN %%%%%%%%%%%%%%%%%%%%%%%%
% rnvar=2*turn_number+3; % change parameters SEED_VALUE
%
% system(['cat newbeam.i | sed -e s/_SEED_VALUE_/' num2str(rnvar) '/ | sed -e s/_pulseE_/' num2str(pulseE) '/ | sed -e s/_angle_/' num2str(angle) '/ > newbeam_tmp.i']);% change paramiters SEED_VALUE in iturn_numbertial CAIN file
% system(['cat ' home_dir '/newbeam.i | sed -e s/_SEED_VALUE_/' num2str(rnvar) '/ | sed -e s/_pulseE_/' num2str(beam_parameters.pulseE) '/ | sed -e s/_bunch_length_initial_/' num2str(beam_parameters.bunch_length_initial) '/ | sed -e s/_number_of_steps_cain_/' num2str(number_of_steps_cain) '/ | sed -e s/_turn_number_/' num2str(turn_number) '/ | sed -e s/_angle_/' num2str(beam_parameters.angle) '/ | sed -e s/_shifting_laser_x_/' num2str(beam_parameters.shifting_laser_x) '/ | sed -e s/_shifting_laser_t_/' num2str(beam_parameters.shifting_laser_t) '/  | sed -e s/_shifting_laser_y_/' num2str(beam_parameters.shifting_laser_y) '/  | sed -e s/_shifting_laser_s_/' num2str(beam_parameters.shifting_laser_s) '/ >  ' DIRECTORY_FOR_CAIN 'newbeam_tmp.i']);% change paramiters SEED_VALUE in iturn_numbertial CAIN file
% system(['cat ' pwd '/newbeam.i | sed -e s/_SEED_VALUE_/' num2str(rnvar) '/ | sed -e s/_pulseE_/' num2str(beam_parameters.pulseE) '/ | sed -e s/_bunch_length_initial_/' num2str(std(beam_property(:,TIME_COORDINATE))) '/ | sed -e s/_number_of_steps_cain_/' num2str(number_of_steps_cain) '/ | sed -e s/_turn_number_/' num2str(turn_number) '/ | sed -e s/_angle_/' num2str(beam_parameters.angle) '/ | sed -e s/_shifting_laser_x_/' num2str(beam_parameters.shifting_laser_x) '/ | sed -e s/_shifting_laser_t_/' num2str(beam_parameters.shifting_laser_t) '/  | sed -e s/_shifting_laser_y_/' num2str(beam_parameters.shifting_laser_y) '/  | sed -e s/_shifting_laser_s_/' num2str(beam_parameters.shifting_laser_s) '/ |     sed -e s/_sigt_/' num2str(beam_parameters.sigt) '/ >  ' DIRECTORY_FOR_CAIN 'newbeam_tmp.i']);% change paramiters SEED_VALUE in iturn_numbertial CAIN file
% home_dir
TIME_COORDINATE=7;%now we start use s(m) for caine
% system(['cat ' home_dir 'newbeam.i | sed -e s/_SEED_VALUE_/' num2str(rnvar) '/ | sed -e s/_pulseE_/' num2str(beam_parameters.pulseE) '/ | sed -e s/_bunch_length_initial_/' num2str(std(beam_property(:,TIME_COORDINATE))) '/ | sed -e s/_number_of_steps_cain_/' num2str(number_of_steps_cain) '/ | sed -e s/_turn_number_/' num2str(turn_number) '/ | sed -e s/_angle_/' num2str(beam_parameters.angle) '/ | sed -e s/_shifting_laser_x_/' num2str(beam_parameters.shifting_laser_x) '/ | sed -e s/_shifting_laser_t_/' num2str(beam_parameters.shifting_laser_t) '/  | sed -e s/_shifting_laser_y_/' num2str(beam_parameters.shifting_laser_y) '/  | sed -e s/_shifting_laser_s_/' num2str(beam_parameters.shifting_laser_s) '/ |     sed -e s/_sigt_/' num2str(beam_parameters.sigt) '/ >  ' DIRECTORY_FOR_CAIN 'newbeam_tmp.i']);% change paramiters SEED_VALUE in iturn_numbertial CAIN file
write_initial_for_cain(turn_number,beam_property)

%system('rm exp.dat'); % delete exp.dat
%%%%%%%%%%%%%%%%%%%%%%%%%% START CAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% system(['cd ' DIRECTORY_FOR_CAIN '; ' home_dir 'cain.ub_illya/cain240/src/cain.exe < newbeam_tmp.i > cain.out']); % start cain
% system(['cd ' DIRECTORY_FOR_CAIN '; ' home_dir 'cain.camilla/cain240/src/cain.exe < newbeam_tmp.i > cain.out']); % start cain
% system(['cd ' DIRECTORY_FOR_CAIN '; ' home_dir 'cain.gcc/cain240/src/cain.exe < newbeam_tmp.i > cain.out']); % start cain
if ismac==1
    system(['cd ' DIRECTORY_FOR_CAIN '; ' home_dir 'cain.mac/cain240/src/cain.exe < newbeam_tmp.i > cain.out']); % start cain
else
    system(['cd ' DIRECTORY_FOR_CAIN '; ' home_dir 'cain.ub_illya/cain240/src/cain.exe < newbeam_tmp.i > cain.out']); % start cain
end

if(DEBUGGING_MODE==1)
    system(['tail -100 ' DIRECTORY_FOR_CAIN 'cain.out']);% it's that see cain
end

beam_property=[];
% beam_property=dlmread([ DIRECTORY_FOR_CAIN 'cain_output_electrons.dat'],'',1,0);%read electrons from cain
%system(['cp -f ' DIRECTORY_FOR_CAIN 'cain_output_photons.dat ' DIRECTORY_FOR_CAIN 'cain_output_photons_' num2str(turn_number) '.dat']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
