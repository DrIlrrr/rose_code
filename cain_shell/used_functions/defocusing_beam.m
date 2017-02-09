function [beam_phasespace] = defocusing_beam(beam_phasespace,defocusing_parameter)% 
% function for defocusing beam on a defocusing factor
% using compensation of longitudinal momentum to save full momentum

  a=defocusing_parameter;%defocusing parameter
 beam_phasespace(:,1)=beam_phasespace(:,1).*a;
 beam_phasespace(:,2)=beam_phasespace(:,2)./a;
 beam_phasespace(:,3)=beam_phasespace(:,3).*a;
 beam_phasespace(:,4)=beam_phasespace(:,4)./a;

  beam_phasespace(:,6)=beam_phasespace(:,6).*sqrt((1+(beam_phasespace(:,2).^2+beam_phasespace(:,4).^2).*(1-1/(a^2))));
