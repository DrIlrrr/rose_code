function write_initial_for_cain(turn_number,beam_property)
global  home_dir DIRECTORY_FOR_CAIN;
global rflags
tilted_name = 'newbeam_tmp';

number_of_macro_electrons=length(beam_property(4,:));
sigz=std(beam_property(4,:));
% angle=rflags.angle*(pi/180); %initial scattered angle [rad]
% pulseE=0.2; %laser puse energy [J]
% sigLr=28/2; % given in [mu m] micro meter like 2 weist w0=28;
% laserwl=515; % laser wavelenth [nm] nano meters
% sigt=1.5; %pulse length [ps]

% shifting_laser_x = 0;  %
% shifting_laser_y = 0;  %
% shifting_laser_s = 0;  %
% shifting_laser_t = 0;%rflags.shifting_laser_t;  %


rnvar=2*turn_number+3;%randi([1, 100]); % change parameters SEED_VALUE

% file filling
filename = sprintf('%s.i',tilted_name);
% fileID = fopen(['CAIN_files/',filename],'w');
fileID = fopen([DIRECTORY_FOR_CAIN filename],'w');
% header
fprintf(fileID,' ALLOCATE  MP=%d;\n',3*number_of_macro_electrons);
fprintf(fileID,' \n');
fprintf(fileID,' SET   photon=1, electron=2, positron=3, mm=1D-3, micron=1D-6, nm=1D-9, mu0=4*Pi*1D-7, \n');
fprintf(fileID,'  psec=1e-12*Cvel,\n');
fprintf(fileID,' \n');
fprintf(fileID,'sigz=%f,\n',sigz);
fprintf(fileID,'ntcut=5,\n');
fprintf(fileID,'laserwl=%f*nm, \n',rflags.laserwl);
fprintf(fileID,'pulseE=%f,\n',rflags.pulseE);
fprintf(fileID,'sigLr=%f*micron,\n',rflags.sigLr);
fprintf(fileID,'w0=2*sigLr,\n');
fprintf(fileID,'rayl=Pi*w0^2/laserwl, \n');
fprintf(fileID,'sigt=%f*psec,\n',rflags.sigt);
fprintf(fileID,'angle=%f,  \n',rflags.angle);
fprintf(fileID,'tdl=1.0,\n');
fprintf(fileID,' powerd=(2*pulseE*Cvel)/[Pi*sigt*Sqrt(2*Pi)*w0^2],\n');
fprintf(fileID,' \n');
fprintf(fileID,' SET MsgLevel=1;\n');
fprintf(fileID,' \n');
fprintf(fileID,' SET Rand=5*%f;\n',rnvar);
fprintf(fileID,' \n');
fprintf(fileID,' BEAM FILE=''exp.dat'';\n');
fprintf(fileID,' \n');
fprintf(fileID,' LASER LEFT, WAVEL=laserwl, POWERD=powerd,\n');
fprintf(fileID,' \n');     
fprintf(fileID,'       TXYS=(%f, %f, %f, %f), \n',rflags.shifting_laser_t,rflags.shifting_laser_x,rflags.shifting_laser_y,rflags.shifting_laser_s);
fprintf(fileID,'       E3=(-Sin(angle),0.0,-Cos(angle)), E1=(0,1,0),\n');
fprintf(fileID,'       RAYLEIGH=(rayl,rayl), SIGT=sigt, GCUTT=ntcut, STOKES=(%f, %f, %f),\n',rflags.STOKES(1),rflags.STOKES(2),rflags.STOKES(3));
fprintf(fileID,'       TDL=(tdl,tdl) ;\n');
fprintf(fileID,' \n');
fprintf(fileID,' \n');
fprintf(fileID,' LASERQED  COMPTON, NPH=0;\n');
fprintf(fileID,' SET MsgLevel=0;  FLAG OFF ECHO;\n');
fprintf(fileID,' SET Smesh=sigt/3;\n');
% fprintf(fileID,' SET emax=1.001*ee, wmax=emax;\n');
fprintf(fileID,' SET  it=0;\n');
fprintf(fileID,' \n');
fprintf(fileID,' PUSH  Time=(-ntcut*(sigt+sigz),ntcut*(sigt+sigz),2000);\n');
fprintf(fileID,'       IF Mod(it,20)=0;\n');
fprintf(fileID,'        PRINT it, FORMAT=(F6.0,''-th time step''); PRINT STAT, SHORT;\n');
fprintf(fileID,'       ENDIF;\n');
fprintf(fileID,'      SET it=it+1;\n');
fprintf(fileID,' ENDPUSH;\n');
fprintf(fileID,' \n');
fprintf(fileID,'DRIFT T=0;\n');
fprintf(fileID,' \n');
fprintf(fileID,'WRITE BEAM, KIND=(electron), FILE=''cain_output_electrons.dat'';\n');
fprintf(fileID,' \n');
fprintf(fileID,'!PRINT STATISTICS, FILE=''ELECRON_STAT.DAT'';\n');
fprintf(fileID,' \n');
fprintf(fileID,'WRITE BEAM, KIND=(photon), FILE=''cain_output_photons_%d.dat'';\n',turn_number);
fprintf(fileID,' \n');
fprintf(fileID,' \n');
fclose(fileID);
%end
save([DIRECTORY_FOR_CAIN 'laser_param.dat'],'rflags')%.laserwl','rflags.pulseE','rflags.sigLr','rflags.sigt',...
%     'rflags.angle','rflags.shifting_laser_t','rflags.shifting_laser_x','rflags.shifting_laser_y','rflags.shifting_laser_s');



