function [el] = global_
% create global flags to run choosen type of simulations
el=struct('NAME_OF_FLAGS',1);

%% electron bunch parameters
el.chargebunch = 250e-12;%Charge per electrons bunch [c] pico->10^-12
el.NUMBER_OF_MACROPARTICLES=1000;%1e5; % Number of macroparticles

el.bunch_length_initial=2.7e-4; % intial bunch length [m]

el.sigma_e_x=17.2*1e-6; % IP vertical electron beam size [m]
el.sigma_e_y=16.4*1e-6; % IP horizontal electron beam size [m]

el.initial_beam_energy_MeV=529.8; % initial energy in [MeV]

el.energy_spread_initial=0.044;%  initial relative energy spread (not in [%])

el.norm_emit_x=0.42*1e-6; %Normilized emittance x [m rad]
el.norm_emit_y=0.44*1e-6; %Normilized emittance y [m rad]
