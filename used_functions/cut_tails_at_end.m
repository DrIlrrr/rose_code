function [beam_1,beam_2]=cut_tails_at_end(beam_1,beam_2,l_max)%,lim_cut)
global ifig gvar

[beam_1_max]=beam_drift(beam_1,l_max);
[beam_2_max]=beam_drift(beam_2,-l_max);

lim_cut=gvar.sigma_cut_at_end*std(beam_1_max(:,5));


% in_10mum_1=find(abs(beam_1_max(:,5))<lim_cut & abs(beam_1_max(:,6))<lim_cut);
in_10mum_1=find((beam_1_max(:,5)).^2+(beam_1_max(:,6)).^2<lim_cut^2);
beam_temp=beam_1(in_10mum_1,:);
beam_1=[];
beam_1=beam_temp;
beam_temp=[];%in_10mum_1=[];


% in_10mum_2=find(abs(beam_2_max(:,5))<lim_cut & abs(beam_2_max(:,6))<lim_cut);
in_10mum_2=find((beam_2_max(:,5)).^2+(beam_2_max(:,6)).^2<lim_cut^2);
beam_temp=beam_2(in_10mum_2,:);
beam_2=[];
beam_2=beam_temp;

beam_temp=[];%in_10mum_2=[];

% close all
% beam_stat('beam_1_cuted_end',beam_1_max(in_10mum_1,:))
% beam_stat('beam_2_cuted_end',beam_2_max(in_10mum_2,:))

% stop
