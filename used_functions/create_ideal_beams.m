function [beam_1,beam_2]=create_ideal_beams
global rflags el

%% create e electron beam
[beam_1] = electron_beam_initial;
% beam_stat(out_folder,beam_1); %save a beam stat


%% create a second beam as the reflection of first one
beam_2=beam_1;
beam_2(:,11)=-beam_2(:,11);
beam_2(:,5)=beam_1(:,6);
beam_2(:,6)=beam_1(:,5);
% [beam_1l]=beam_drift(beam_1,0);
%

%% put beams around z=0
beam_1(:,4)=beam_1(:,4)-max(beam_1(:,4));
beam_2(:,4)=beam_2(:,4)-min(beam_2(:,4));

% beam_1(:,4)=beam_1(:,4).*5;
% beam_2(:,4)=beam_2(:,4)./5;






