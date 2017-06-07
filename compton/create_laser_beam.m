function [beam_phasespace] = create_laser_beam
global el
%% Create a laser beam
% LASER parameteres
NUMBER_OF_MACROPARTICLES=el.NUMBER_OF_MACROPARTICLES;
angle_deg=8;% in degree
angle=angle_deg*(pi/180); %initial scattered angle [rad]
pulseE=0.4; %laser puse energy [J]
sigLr=14e-6; % given in [m] micro meter like 2 weist w0=28;
laserwl=515e-9; % laser wavelenth [m] nano meters
sigt=1.5e-12; %pulse length [s]2.7e-4/3e8;%
shifting_laser_x = 0;  %
shifting_laser_y = 0;  %
shifting_laser_s = 0;  %
shifting_laser_t = 0;%shifting_laser_t;
STOKES=[0 0 1];% linear

tr_per=1e-3;
%%
Xdeviation =sigLr*randn(1,NUMBER_OF_MACROPARTICLES);% vector of initial vertical laser beam size m
Ydeviation =sigLr*randn(1,NUMBER_OF_MACROPARTICLES);% vector of initial horizontal laser beam size m

SPEED_OF_LIGHT=3e8;
h=2*pi*1.054e-34;%Planc const [J*s]
h_ev=4.135e-15;%Planc const [eV*s]
initial_laser_energy=h_ev*SPEED_OF_LIGHT/laserwl;
delta_nu=1/(2*pi*sigt);
laser_BW=h_ev*delta_nu;
randn_for_all=randn(1,NUMBER_OF_MACROPARTICLES);
randn_for_all=randn_for_all-mean(randn_for_all);

energy_deviation= initial_laser_energy.*(1+laser_BW.*randn_for_all);% vector of initial electron beam energy
bunch_length = sigt*SPEED_OF_LIGHT*randn(1,NUMBER_OF_MACROPARTICLES);%-mean(bunch_length_NO_OFFSET); % vector of initial bunch length [m]


photons_number=pulseE/((h*SPEED_OF_LIGHT)/laserwl);
weight = photons_number/NUMBER_OF_MACROPARTICLES;%Number of real particles in one macroparticle
weight0 = ones(1,NUMBER_OF_MACROPARTICLES)*weight;

momentum=sqrt(energy_deviation.^2);%/SPEED_OF_LIGHT;normalezd for cain % Momentum full aka Pi [Mev/c]
Pzi = momentum./sqrt(1+(tr_per.*randn_for_all).^2);% Longitudinal momentum[Mev/c]
Pti=sqrt(momentum.^2-Pzi.^2);
momangle = 2.*pi.*rand(1,NUMBER_OF_MACROPARTICLES);
Pxi = Pti.*cos(momangle); % Momentum projections[Mev/c]
Pyi = Pti.*sin(momangle);% Momentum projections[Mev/c]


% Creating polarizations
Sx = STOKES(1).*ones(1,NUMBER_OF_MACROPARTICLES);% X polarizations
Sy = STOKES(2).*ones(1,NUMBER_OF_MACROPARTICLES);% Y polarizations
Ss = STOKES(3).*ones(1,NUMBER_OF_MACROPARTICLES);% s polarizations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%data for 2 first columns in cain standart data file
K=2;% for cain data output
genname=1;% for cain data output
K0 = ones(1,NUMBER_OF_MACROPARTICLES)*K;
genname0 = ones(1,NUMBER_OF_MACROPARTICLES)*genname;
%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss
beam_phasespace=[K0;genname0;weight0;bunch_length;Xdeviation;Ydeviation;0*bunch_length;energy_deviation;Pxi;Pyi;Pzi;Sx;Sy;Ss;]';
