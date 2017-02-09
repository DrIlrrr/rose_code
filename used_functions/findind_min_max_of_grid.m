function [x_max,x_min,y_max,y_min,z_max,z_min,l_max]=findind_min_max_of_grid(beam_1,beam_2)
% this fuction find a min a max of XYZ for two beams what stay around z=0
%% create a greed

%% traditional
% z_min=min([min(beam_1(:,4)) min(beam_2(:,4))]);
% z_max=max([max(beam_1(:,4)) max(beam_2(:,4))]);
% % %l_max distance what take the farest particle to cross z=0
% % l_max=max([abs(max(beam_1(:,4))) abs(max(beam_2(:,4))) abs(min(beam_1(:,4))) abs(min(beam_2(:,4)))]);
%
% l_max=(abs(z_max)+abs(z_min))/2;
%
% [beam_1_max]=beam_drift(beam_1,l_max);
% [beam_2_max]=beam_drift(beam_2,-l_max);
% x_min=min([min(beam_1_max(:,5)) min(beam_2_max(:,5))]);
% x_max=max([max(beam_1_max(:,5)) max(beam_2_max(:,5))]);
% y_min=min([min(beam_1_max(:,6)) min(beam_2_max(:,6))]);
% y_max=max([max(beam_1_max(:,6)) max(beam_2_max(:,6))]);
% % x_min=min([min(beam_1_max(:,5)) min(beam_2_max(:,5)) min(beam_1(:,5)) min(beam_2(:,5))]);
% % x_max=max([max(beam_1_max(:,5)) max(beam_2_max(:,5)) max(beam_1(:,5)) max(beam_2(:,5))]);
% % y_min=min([min(beam_1_max(:,6)) min(beam_2_max(:,6)) min(beam_1(:,6)) min(beam_2(:,6))]);
% % y_max=max([max(beam_1_max(:,6)) max(beam_2_max(:,6)) max(beam_1(:,6)) max(beam_2(:,6))]);
% %
% %
% %
% % z_min=min([min(beam_1_max(:,4)) min(beam_2_max(:,4)) min(beam_1(:,4)) min(beam_2(:,4))]);
% % z_max=max([max(beam_1_max(:,4)) max(beam_2_max(:,4)) max(beam_1(:,4)) max(beam_2(:,4))]);

%% just for one beam
% must be faster for simulations with initial alpha
z_min=min([min(beam_1(:,4)) min(beam_2(:,4))]);
z_max=max([max(beam_1(:,4)) max(beam_2(:,4))]);
% %l_max distance what take the farest particle to cross z=0
% l_max=max([abs(max(beam_1(:,4))) abs(max(beam_2(:,4))) abs(min(beam_1(:,4))) abs(min(beam_2(:,4)))]);

l_max=(abs(z_max)+abs(z_min))/2;

[beam_1_max]=beam_drift(beam_1,l_max);
[beam_2_max]=beam_drift(beam_2,-l_max);
x_min= min(beam_1_max(:,5));
x_max= max(beam_1_max(:,5));
y_min= min(beam_1_max(:,6));
y_max= max(beam_1_max(:,6));
