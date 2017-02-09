%%%%%%%%%%%%%%%%DO NOT CHANGE ANYTHING BELOW (AFTER THIS LINE)%%%%%%%%%%%
% Creating the electron beam 
%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss
 function[beam_phasespace] = FEL_gamma_initial
 global rflag
% global BASE_DIRECTORY beam_parameters home_dir;
% clear all; close all; clc;
BASE_DIRECTORY=[pwd]
for_plot=0;
ifig=1;

SPEED_OF_LIGHT=3e8; % velocity of light [m/s]
bunch_length_initial=2.7e-4; % intial bunch length [m]

%% not change

 Emass = 0;% photon energy of rest
% gamma=el.initial_beam_energy_MeV*1e6/Emass;
% 
% emit_x=el.norm_emit_x/sqrt(gamma^2-1);
% emit_y=el.norm_emit_y/sqrt(gamma^2-1);
% 
% el.Betax=el.sigma_e_x^2/emit_x;
% el.Betay=el.sigma_e_y^2/emit_y;

initial_beam_energy=1e4; %[eV]
energy_spread_initial=0*1e-3;
sigma_p_x=1e-6;
sigma_p_y=1e-6;
Betax=0.1;
Betay=0.1;

NUMBER_OF_MACROPARTICLES=1e4; % Number of macroparticles
number_photons = 1e13;% Number photons in bunch
weight = number_photons/NUMBER_OF_MACROPARTICLES;%Number of real particles in one macroparticle

randn_for_all=randn(1,NUMBER_OF_MACROPARTICLES);
randn_for_all=randn_for_all-mean(randn_for_all);

energy_deviation= initial_beam_energy.*(1+energy_spread_initial.*randn_for_all);% vector of initial electron beam energy
bunch_length = bunch_length_initial*randn(1,NUMBER_OF_MACROPARTICLES);%-mean(bunch_length_NO_OFFSET); % vector of initial bunch length [m]

Xdeviation =sigma_p_x*randn(1,NUMBER_OF_MACROPARTICLES);% vector of initial vertical electron beam size m
Ydeviation =sigma_p_y*randn(1,NUMBER_OF_MACROPARTICLES);% vector of initial horizontal electron beam size m

%data for 2 first columns in cain standart data file
K=2;% for cain data output
genname=1;% for cain data output
K0 = ones(1,NUMBER_OF_MACROPARTICLES)*K;
genname0 = ones(1,NUMBER_OF_MACROPARTICLES)*genname;
weight0 = ones(1,NUMBER_OF_MACROPARTICLES)*weight;
t =zeros(1,NUMBER_OF_MACROPARTICLES);% for time of work cain
% Creating momentum
momentum = sqrt(energy_deviation.^2-Emass^2);%/SPEED_OF_LIGHT;normalezd for cain % Momentum full aka Pi [Mev/c]

%  momangle = rand(1,NUMBER_OF_MACROPARTICLES);
X_prime=(sigma_p_x/Betax)*randn(1,NUMBER_OF_MACROPARTICLES);%.*cos(momangle);
Y_prime=(sigma_p_y/Betay)*randn(1,NUMBER_OF_MACROPARTICLES);%.*sin(momangle);

Pzi = momentum./sqrt(1+X_prime.^2+Y_prime.^2);% Longitudinal momentum[Mev/c]
Pxi = X_prime.*Pzi; % Momentum projections[Mev/c]
Pyi = Y_prime.*Pzi;% Momentum projections[Mev/c]

% Creating polarizations
Sx = zeros(1,NUMBER_OF_MACROPARTICLES);% X polarizations
Sy = zeros(1,NUMBER_OF_MACROPARTICLES);% Y polarizations
Ss = zeros(1,NUMBER_OF_MACROPARTICLES);% s polarizations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  beam_property_initial=transpose([K0;genname0;weight0;0*bunch_length;Xdeviation;Ydeviation;bunch_length;energy_deviation;Pxi;Pyi;Pzi;Sx;Sy;Ss;]);
%   colons=[K0;genname0;weight0;bunch_length;Xdeviation;Ydeviation;0*bunch_length;energy_deviation;Pxi;Pyi;Pzi;Sx;Sy;Ss];
%   fid = fopen('exp.dat','w');
%   fprintf(fid,' %i    %i       %1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e \n',colons);
%   fclose(fid);
%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss
beam_phasespace=[K0;genname0;weight0;bunch_length;Xdeviation;Ydeviation;0*bunch_length;energy_deviation;Pxi;Pyi;Pzi;Sx;Sy;Ss;]';  