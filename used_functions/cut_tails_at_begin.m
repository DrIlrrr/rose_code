function [beam_1,beam_2]=cut_tails_at_begin(beam_1,beam_2,lim_cut)
global ifig


in_10mum_1=find((beam_1(:,5)).^2+(beam_1(:,6)).^2<lim_cut^2);
beam_temp=beam_1(in_10mum_1,:);
beam_1=[];
beam_1=beam_temp;
beam_temp=[];%in_10mum_1=[];


% in_10mum_2=find(abs(beam_2_max(:,5))<lim_cut & abs(beam_2_max(:,6))<lim_cut);
in_10mum_2=find((beam_2(:,5)).^2+(beam_2(:,6)).^2<lim_cut^2);
beam_temp=beam_2(in_10mum_2,:);
beam_2=[];
beam_2=beam_temp;

beam_temp=[];%in_10mum_2=[];

% close all 
% beam_stat('beam_1_cuted_end',beam_1_max(in_10mum_1,:))
% beam_stat('beam_2_cuted_end',beam_2_max(in_10mum_2,:))

% stop