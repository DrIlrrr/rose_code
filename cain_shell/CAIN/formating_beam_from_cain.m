function [beam_phasespace] = formating_beam_from_cain(beam_property)
global  beam_parameters 
%%%%%%%%%%%%%%%%%  formating  beam data from cain %%%%%%%%%%%%%%%%%%%%%%%%%%
WEIGHT=3;
TIME_COORDINATE=7;%now we start use s(m) for caine
X_COORDINATE=5;
Y_COORDINATE=6;
ENERGY_OF_PARTICLE=8;
X_MOMENTUM=9;
Y_MOMENTUM=10;
S_MOMENTUM=11;
% POLARISATION: 12 13 14
N_COMPTON_HIT=15;
TURN_LAST_COMPTON_HIT=16;

beam_phasespace(1,:)=(beam_property(:,X_COORDINATE))';

beam_phasespace(2,:)=(beam_property(:,X_MOMENTUM)./(beam_property(:,S_MOMENTUM)))';

beam_phasespace(3,:)=(beam_property(:,Y_COORDINATE))';

beam_phasespace(4,:)=(beam_property(:,Y_MOMENTUM)./(beam_property(:,S_MOMENTUM)))';

beam_phasespace(5,:)=beam_property(:,TIME_COORDINATE)';
% if(beam_parameters.no_S_propagaionts_in_CAIN==1)
%
%  beam_phasespace(5,:)=beam_phasespace(5,:)-mean(beam_phasespace(5,:));
% end


beam_phasespace(6,:)=((beam_property(:,ENERGY_OF_PARTICLE)-beam_parameters.initial_beam_energy)/beam_parameters.initial_beam_energy)';


beam_phasespace(7,:)=(beam_property(:,1))';
beam_phasespace(8,:)=(beam_property(:,2))';
beam_phasespace(9,:)=(beam_property(:,WEIGHT))';
beam_phasespace(10,:)=(beam_property(:,4))';

beam_phasespace(11,:)=(beam_property(:,12))';
beam_phasespace(12,:)=(beam_property(:,13))';
beam_phasespace(13,:)=(beam_property(:,14))';
%     beam_phasespace(14,:)=(beam_property(:,15))';
%     beam_phasespace(15,:)=(beam_property(:,16))';

%number_of_scattered_electrons=beam_property(:,2)-1;

% number_of_scattered_electrons=beam_phasespace(8,:)-1;
% 
% beam_phasespace(8,:)=1;

% indices=find(number_of_scattered_electrons>0);
% beam_phasespace(14,indices)=beam_phasespace(14,indices)+1;
% beam_phasespace(15,indices)=turn_number;


