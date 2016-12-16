function [turned_beam]=turn_beam(alpha,beam)

%turn a beam on angle alpha around x-axes 
%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss

%% turn around X-axes
turn_matrix=...
    [1 0 0;...
    0 cos(alpha) -sin(alpha);...
    0 sin(alpha) cos(alpha)]
%% turn around Y-axes
%  turn_matrix=...
%      [cos(alpha) 0 sin(alpha);...
%      0 1 0;...
%      -sin(alpha) 0 cos(alpha)]
%% turn around Z-axes
%  turn_matrix=...
%      [cos(alpha) -sin(alpha) 0;...
%       sin(alpha) cos(alpha) 0;...
%      0 0 1]
%%
gvec=[beam(:,5)'; beam(:,6)'; beam(:,4)'];

mom_vec=[beam(:,9)'; beam(:,10)'; beam(:,11)'];


beam_temp=turn_matrix*gvec(1:3,:);
beam_temp_m=turn_matrix*mom_vec(1:3,:);
turned_beam=beam;

turned_beam(:,5)=beam_temp(1,:);
turned_beam(:,6)=beam_temp(2,:);
turned_beam(:,4)=beam_temp(3,:);
turned_beam(:,9)=beam_temp_m(1,:);
turned_beam(:,10)=beam_temp_m(2,:);
turned_beam(:,11)=beam_temp_m(3,:);
