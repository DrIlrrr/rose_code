clear all; close all; clc;

ev=load('rs_grid.dat','ascii')';
cs=load('w_grid.dat','ascii')';
cross=load('cs_grid.dat','ascii')';

figure(1)
mesh((cs(:,1)),ev(1,:),cross')
xlabel('\theta')
ylabel('E_{cm}')
view(44,30)