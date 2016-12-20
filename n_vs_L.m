clear all; close all; clc;

x_1=[0 1 2 3 4 5 10];
y_11=[3.44e7 3.09e7 2.42e7 1.80e7 1.61e7 1.20e7 0];
y_51=[3.66e7 3.25e7 2.56e7 1.92e7 1.66e7 1.27e7 6.56e6];

% L_TPP=[0 1e-3 2e-3 4e-3 5e-3 6e-3 8e-3 1e-2];
% N_TPP=[5.55e-5 4.68e-5 3.87e-5 2.7e-5 2.34e-5 1.99e-5 1.48e-5 1.12e-5];



figure(1)
set(gca,'FontSize',16)
hold on
plot(x_1,y_11,'-*','LineWidth',2)
plot(x_1,y_51,'-*','LineWidth',2)

hold off
grid on
legend('BW','TPP','FontSize',20)
xlabel('L mm','FontSize',20)
ylabel('N events','FontSize',20)
fname = ['num_vs_L.png'];
print('-r300','-dpng', fname);
