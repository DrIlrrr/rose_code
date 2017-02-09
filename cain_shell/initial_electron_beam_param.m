function [el] = initial_electron_beam_param
% create global flags to run choosen type of simulations
el=struct('NAME_OF_FLAGS',1);

%% electron bunch parameters
el.chargebunch = 40e-12;%Charge per electrons bunch [c] pico->10^-12
el.NUMBER_OF_MACROPARTICLES=1e5; % Number of macroparticles

el.bunch_length_initial=3e-5; % intial bunch length [m]

el.sigma_e_x=7.7*1e-6; % IP vertical electron beam size [m]
el.sigma_e_y=7.7*1e-6; % IP horizontal electron beam size [m]

el.initial_beam_energy_MeV=7e3; % initial energy in [MeV]

el.energy_spread_initial=0*2e-4;%  initial relative energy spread (not in [%])

el.norm_emit_x=0.08*1e-6; %Normilized emittance x [m rad]
el.norm_emit_y=0.08*1e-6; %Normilized emittance y [m rad]
