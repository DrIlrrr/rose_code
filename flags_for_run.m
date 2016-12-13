
function [rflags] = flags_for_run
% create global flags to run choosen type of simulations
rflags=struct('NAME_OF_FLAGS',1);

% rflags.NAME_OF_PARAMS = ['default']; %Number of turns
YES=1; yes=1; NO=0; no=0;
rflags.name_of_job=[''];
rflags.PLOTS =yes;

rflags.e_gamma=0;%for muon prodaction
rflags.gamma_gamma=0;
rflags.breit_wheeler=0;
rflags.compton=0;%Compton Back Scattering
rflags.TPP=0;%triplet pair production
rflags.moller=0;


