function write_initial_for_cain(turn_number,beam_property)
global  home_dir DIRECTORY_FOR_CAIN;
global rflags
tilted_name = 'newbeam_tmp';


angle=7.5*pi/180; %initial scattered angle [rad]
pulseE=0.2; %laser puse energy [J]
sigt=3.5; %pulse length [ps]
% sigz=std(beam_property(4,:));
shifting_laser_x = 0;  %
shifting_laser_y = 0;  %
shifting_laser_s = 0;  %
shifting_laser_t = rflags.shifting_laser_t;  %


rnvar=2*turn_number+3;%randi([1, 100]); % change parameters SEED_VALUE

% file filling
filename = sprintf('%s.i',tilted_name);
% fileID = fopen(['CAIN_files/',filename],'w');
fileID = fopen([DIRECTORY_FOR_CAIN filename],'w');
% header
fprintf(fileID,' ALLOCATE  MP=2500000;\n');
fprintf(fileID,' FLAG OFF ECHO;\n');
fprintf(fileID,' HEADER  ''Laser-Compton;'';\n');
fprintf(fileID,' SET   photon=1, electron=2, positron=3,\n');
fprintf(fileID,'mm=1D-3, micron=1D-6, nm=1D-9, mu0=4*Pi*1D-7, psec=1e-12*Cvel,;\n');
fprintf(fileID,'\n');
fprintf(fileID,'SET sigz=1*psec;\n');
fprintf(fileID,'! Laser parameters\n');
fprintf(fileID,' SET wl=515*nm, pulseE=%f, wlbar=wl/(2*Pi),\n',pulseE);
% fprintf(fileID,'! LASER pulse energy 1 J \n');
fprintf(fileID,'   eph=Hbarc/wlbar,   !  laser photon energy in eV\n');
fprintf(fileID,'   sigLr=30*micron/2.0, rayl=2/wlbar*sigLr^2, w0=2*sigLr, ! Waist is 2sigLr  \n');
fprintf(fileID,'   sigt=%f*psec,\n',sigt);
fprintf(fileID,'   angle=%f , tdl=1.0,\n',angle);
fprintf(fileID,'   p0=(2*pulseE*Cvel) / [Pi*sigt*Sqrt(2*Pi)*w0^2],\n');
fprintf(fileID,'   xisq=p0*mu0*Cvel*(wlbar/Emass)^2,   xi=Sqrt(xisq);\n');
fprintf(fileID,'\n');
fprintf(fileID,'PRINT xi, FORMAT=("Maximum Xi parameter =",F10.5);\n');
fprintf(fileID,'PRINT eph, FORMAT=("Laser photon energy =",F10.5," eV");\n');
fprintf(fileID,'PRINT rayl, FORMAT=("ZRai1l =",F10.5);\n');
fprintf(fileID,'PRINT rayl2, FORMAT=("ZRail2 =",F10.5);\n');
fprintf(fileID,'PRINT p0, FORMAT=("PowerDensity = ",D10.3," W/m2");\n');
fprintf(fileID,'\n');
fprintf(fileID,' !SET MsgLevel=1;\n');
fprintf(fileID,' SET Rand= %f;\n',rnvar);
fprintf(fileID,'\n');
fprintf(fileID,'! Define initial electron beam at the focal point\n');
fprintf(fileID,' BEAM FILE=''exp.dat'';\n');
fprintf(fileID,'\n');
fprintf(fileID,'! Define LaserQED parameters\n');
fprintf(fileID,' LASERQED  COMPTON, CIRCULARPOL, NPH=0; \n');
fprintf(fileID,'\n');
fprintf(fileID,' SET MsgLevel=0;\n');
fprintf(fileID,'\n');
fprintf(fileID,' LASER  LEFT, WAVELENGTH=wl, POWERDENSITY=p0, \n');
fprintf(fileID,' TXYS=(%f, %f, %f, %f), \n',shifting_laser_t,shifting_laser_x,shifting_laser_y,shifting_laser_s);
fprintf(fileID,'\n');
fprintf(fileID,' E3=(-Sin(angle),0,-Cos(angle)), E1=(0,1,0), STOKES=(0,1,0),\n');
fprintf(fileID,'\n');
fprintf(fileID,' SIGT=sigt,\n');
fprintf(fileID,' RAYLEIGH=(rayl,rayl),GCUTT=5;\n');
fprintf(fileID,' SET  it=0;\n');
fprintf(fileID,' PUSH   Time=(-5*sigz-5*sigt, 5*sigz+5*sigt, 200);\n');
fprintf(fileID,'   IF Mod(it,20)=0;\n');
fprintf(fileID,'        PRINT it, FORMAT=(F6.0,''-th time step''); PRINT STAT, SHORT;\n');
fprintf(fileID,'      ENDIF;\n');
fprintf(fileID,'      SET it=it+1;\n');
fprintf(fileID,'   ENDPUSH;\n');
fprintf(fileID,'   CLEAR LASER;\n');
fprintf(fileID,'   PRINT (1), FORMAT=("==== After",I2,"-th Interaction ====");\n');
fprintf(fileID,'   PRINT STAT;\n');
fprintf(fileID,'\n');
fprintf(fileID,'!  Pull all particles to the IP\n');
fprintf(fileID,' DRIFT S=0;\n');
fprintf(fileID,'\n');
fprintf(fileID,' ! Write particle data into a file s\n');
fprintf(fileID,'WRITE BEAM, KIND=(electron), FILE=''cain_output_electrons.dat'';\n');
fprintf(fileID,'\n');
fprintf(fileID,'PRINT STATISTICS, FILE=''ELECRON_STAT.DAT'';\n');
fprintf(fileID,'\n');
fprintf(fileID,'WRITE BEAM, KIND=(photon), FILE=''cain_output_photons.dat'';\n');
fprintf(fileID,'  \n');
fprintf(fileID,' STOP;\n');
fprintf(fileID,' END;\n');
fprintf(fileID,'\n');
fprintf(fileID,'\n');
fprintf(fileID,'\n');


fclose(fileID);
%end


