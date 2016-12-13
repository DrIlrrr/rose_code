function [beam_1]=beam_drift(beam_1,l)
%  propagate beam on distance l [m]
% x_1=beam_1(:,5);
% y_1=beam_1(:,6);
% z_1=beam_1(:,4);
% px_1=beam_1(:,9);
% py_1=beam_1(:,10);
% pz_1=beam_1(:,11);
% l=2e-3;
% x_l=x_1+(px_1./pz_1).*l;
% y_l=y_1+(py_1./pz_1).*l;
% z_l=z_1+l;


beam_1(:,1:3)=beam_1(:,1:3);
beam_1(:,4)=beam_1(:,4)+l;
beam_1(:,5)=beam_1(:,5)+(beam_1(:,9)./beam_1(:,11)).*l;
beam_1(:,6)=beam_1(:,6)+(beam_1(:,10)./beam_1(:,11)).*l;
beam_1(:,7:14)=beam_1(:,7:14);

