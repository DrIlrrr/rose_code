function [b_emit] = beam_emit(beam_property)

WEIGHT=3;
TIME_COORDINATE=4;%now we start use s(m) for caine
X_COORDINATE=5;
Y_COORDINATE=6;
ENERGY_OF_PARTICLE=8;
X_MOMENTUM=9;
Y_MOMENTUM=10;
S_MOMENTUM=11;
% POLARISATION: 12 13 14
N_COMPTON_HIT=15;
TURN_LAST_COMPTON_HIT=16;



phasespace(1,:)=(beam_property(X_COORDINATE,:))';

phasespace(2,:)=(beam_property(X_MOMENTUM,:)./(beam_property(S_MOMENTUM,:)))';

phasespace(3,:)=(beam_property(Y_COORDINATE,:))';

phasespace(4,:)=(beam_property(Y_MOMENTUM,:)./(beam_property(S_MOMENTUM,:)))';

phasespace(5,:)=beam_property(TIME_COORDINATE,:)';

phasespace(6,:)=((beam_property(ENERGY_OF_PARTICLE,:)-mean(beam_property(ENERGY_OF_PARTICLE,:)))/mean(beam_property(ENERGY_OF_PARTICLE,:)))';




% std
sig =cov(phasespace(:,:)'); % Sigma matrix 6*6
b_sig =sqrt(diag(sig));     % std of the 6 variables


% emittances
i=1;emitx=sqrt(det(sig(i:i+1,i:i+1)));
i=3;emity=sqrt(det(sig(i:i+1,i:i+1)));
i=5;emits=sqrt(det(sig(i:i+1,i:i+1)));
b_emit=[emitx ; emity ; emits];% non normalized [m rad] emittance of the 3 subspaces 


fprintf('Emittances:\n')
fprintf('Emit X = %10.5e\n',emitx)
fprintf('Emit Y = %10.5e\n',emity)
fprintf('Emit S = %10.5e\n',emits)




