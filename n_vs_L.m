clear all; close all; clc;

L_BW=[0 0.1e-3 0.5e-3 0.8e-3 0.9e-3 1e-3 1.1e-3 1.2e-3 1.3e-3 1.5e-3 2e-3 3e-3 4e-3 5e-3 6e-3 8e-3 1e-2];
N_BW=[7.02e-4 6.67e-4 5.59e-4 4.78e-4 4.55e-4 4.37e-4 4.17e-4 4.10e-4 3.90e-4 3.61e-4 3.02e-4 2.10e-4 1.5e-4 1.16e-4 9.10e-5 6.1e-5 4.28e-5]; 

L_TPP=[0 1e-3 2e-3 4e-3 5e-3 6e-3 8e-3 1e-2];
N_TPP=[5.55e-5 4.68e-5 3.87e-5 2.7e-5 2.34e-5 1.99e-5 1.48e-5 1.12e-5];



figure(1)
set(gca,'FontSize',16)
hold on
plot(L_BW*1e3,N_BW,'-*','LineWidth',2)
plot(L_TPP*1e3,N_TPP,'-*','LineWidth',2)

hold off
grid on
legend('BW','TPP','FontSize',20)
xlabel('L mm','FontSize',20)
ylabel('N events','FontSize',20)
fname = ['num_vs_L.png'];
print('-r300','-dpng', fname); 
