!  ThomX Anneau Compton interaction
 HEADER  'ThomX_Anneau';
 ALLOCATE MP=1000000;
 SET   photon=1, electron=2, positron=3, mm=1D-3, micron=1D-6, nm=1D-9, mu0=4*Pi*1D-7, psec=1e-12*Cvel,
   ee=50D6,  
gamma=ee/Emass,  
an=1e-9/Echarge,  
emitx=4.9D-6/gamma, 
emity=4.9D-6/gamma,
sigz=20*psec,
sige=
0.6D-2
, 
betax=100*mm,
betay=100*mm,
!sigx=Sqrt(emitx*betax), 
!betay=Sqrt(emity*betay),
ntcut=5.0,
laserwl=1.06*micron, 
pulseE=30D-3, 
!omegal=Hbarc/lambar, 
sigLr=40*micron,
w0=2*sigLr,
rayl=Pi*w0^2/laserwl, 
sigt=1*psec,
angle=
0
,  
tdl=1.0,
powerd= (pulseE*Cvel) / [Pi*sigt*Sqrt(2*Pi)*w0^2],
!  	xisq=powerd*mu0*Cvel*(lambar/Emass)^2,   
!	xi=Sqrt(xisq),
!   	lambda=4*omegal*ee/Emass^2;

 SET MsgLevel=1;

 BEAM  RIGHT, KIND=2, NP=600000, AN=an, E0=ee,
   TXYS=(0,0,0,0),  GCUTT=ntcut, SIGE=sige,
   BETA=(betax,betay), EMIT=(emitx,emity), SIGT=sigz, SPIN=(0,0,0);

 LASER LEFT, WAVEL=laserwl, POWERD=powerd,
      TXYS=(0,0,0,0),
      E3=(-Sin(angle),0.0,-Cos(angle)), E1=(1,0,0), 
      RAYLEIGH=(rayl,rayl), SIGT=sigt, GCUTT=ntcut, STOKES=(0,1,0),
      TDL=(tdl,tdl) ;

! LASERQED  COMPTON, LINEARPOL, NPH=5, XIMAX=1.1*xi, LAMBDAMAX=1.1*lambda,
!   PMAX=0.5;


 LASERQED  COMPTON, NPH=0;
 SET MsgLevel=0;  FLAG OFF ECHO;
 SET Smesh=sigt/3;
 SET emax=1.001*ee, wmax=emax;
 SET  it=0;
 PRINT CPUTIME;
 PUSH  Time=(-ntcut*(sigt+sigz),ntcut*(sigt+sigz),200);
      IF Mod(it,20)=0;
        PRINT it, FORMAT=(F6.0,'-th time step'); PRINT STAT, SHORT;
      ENDIF;
      SET it=it+1;
 ENDPUSH;
 PRINT CPUTIME;
!  Pull all particles to the IP
 DRIFT S=0;
!  Write particle data onto a file for the job for IP
! WRITE BEAM, KIND=(electron,photon), FILE='temp.dat';
!  Store variables
 STORE FILE='temp2.dat';
 PRINT STAT;

WRITE BEAM, KIND=(electron), FILE='example_eletrons.dat';

PRINT STATISTICS, FILE='ELECRON_STAT.DAT';

! STOP;
