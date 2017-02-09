function beam_stat(nome_dir,beam_in)
% Make plot and save beam stat
global ifig rflags save_dir_start
%if rflags.PLOTS ==1;
save_dir_plot = [save_dir_start nome_dir ];
mkdir(save_dir_plot);


%beam_in=load(['exp.dat']);can be used for reading "CAIN" output

%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss
WEIGHT=3;
TIME_COORDINATE=4;%now we start use s(m) for caine
X_COORDINATE=5;
Y_COORDINATE=6;
ENERGY_OF_PARTICLE=8;
X_MOMENTUM=9;
Y_MOMENTUM=10;
S_MOMENTUM=11;
% % POLARISATION: 12 13 14
% N_COMPTON_HIT=15;
% TURN_LAST_COMPTON_HIT=16;

electron_angle=sign(beam_in(:,9)).*atan(sqrt(beam_in(:,9).^2+beam_in(:,10).^2)./beam_in(:,11));

ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
subplot 221
set(gca,'FontSize',12)
hist(1e6*beam_in(:,X_COORDINATE),50)
grid on
xlabel('x [\mum]')
title({['STD = ' num2str(std(1e6*beam_in(:,X_COORDINATE))) ' \mum'];['MEAN = ' num2str(mean(1e6*beam_in(:,X_COORDINATE))) ' \mum']},'FontSize',9)
subplot 222
set(gca,'FontSize',11)
hist(1e6*beam_in(:,Y_COORDINATE),50)
grid on
xlabel('y [\mum]')
title({['STD = ' num2str(std(1e6*beam_in(:,Y_COORDINATE))) ' \mum'];['MEAN = ' num2str(mean(1e6*beam_in(:,Y_COORDINATE))) ' \mum']},'FontSize',9)
subplot 223
set(gca,'FontSize',11)
hist(1e3*beam_in(:,TIME_COORDINATE),50)
grid on
xlabel('S [mm]')
title(['STD = ' num2str(std(1e3*beam_in(:,TIME_COORDINATE))) ' mm'],'FontSize',9)
subplot 224
set(gca,'FontSize',11)
hist(beam_in(:,ENERGY_OF_PARTICLE)/1e6,50)
grid on
xlabel('E_{tot} [MeV]')
title(['STD = ' num2str(std(beam_in(:,ENERGY_OF_PARTICLE)/1e6)) ' MeV ' ' Mean = ' num2str(mean(beam_in(:,ENERGY_OF_PARTICLE)/1e6))],'FontSize',9)
filename = [save_dir_plot '/EXP_fig_' num2str(ifig)];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
%print('-depsc', fname);
print('-r300','-dpng', fname2);


if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
ifig=ifig+1;

subplot 221
set(gca,'FontSize',12)
hist(1e-6*beam_in(:,X_MOMENTUM),50)
grid on
xlabel('Px [MeV]')
title({['STD = ' num2str(std(1e-6*beam_in(:,X_MOMENTUM))) ' MeV'];['MEAN = ' num2str(mean(1e-6*beam_in(:,X_MOMENTUM))) ' MeV']},'FontSize',9)
subplot 222
set(gca,'FontSize',11)
hist(1e-6*beam_in(:,Y_MOMENTUM),50)
grid on
xlabel('Py [MeV]')
title({['STD = ' num2str(std(1e-6*beam_in(:,Y_MOMENTUM))) ' MeV'];['MEAN = ' num2str(mean(1e-6*beam_in(:,Y_MOMENTUM))) ' MeV']},'FontSize',9)
subplot 223
set(gca,'FontSize',11)
hist(1e-6*beam_in(:,S_MOMENTUM),50)
grid on
xlabel('Pz [MeV]')
title(['STD = ' num2str(std(1e-6*beam_in(:,S_MOMENTUM))) ' MeV' ' Mean = ' num2str(mean(beam_in(:,S_MOMENTUM)/1e6))],'FontSize',9)
subplot 224
set(gca,'FontSize',11)
hist(beam_in(:,ENERGY_OF_PARTICLE)/1e6,50)
grid on
xlabel('E_{tot} [MeV]')
title(['STD = ' num2str(std(beam_in(:,ENERGY_OF_PARTICLE)/1e6)) ' MeV ' ' Mean = ' num2str(mean(beam_in(:,ENERGY_OF_PARTICLE)/1e6))],'FontSize',9)
filename = [save_dir_plot '/EXP_fig_' num2str(ifig)];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
%print('-depsc', fname);
print('-r300','-dpng', fname2);

if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
ifig=ifig+1;
set(gca,'FontSize',16)
plot(1e3*beam_in(:,TIME_COORDINATE),beam_in(:,ENERGY_OF_PARTICLE)/1e6,'.r')
grid on
ylabel('E [MeV]')
xlabel('S [mm]')
title({['STD = ' num2str(std(1e3*beam_in(:,TIME_COORDINATE))) ' mm'];['STD = ' num2str(std(beam_in(:,ENERGY_OF_PARTICLE)/1e6)) ' MeV ' ' Mean = ' num2str(mean(beam_in(:,ENERGY_OF_PARTICLE)/1e6))]},'FontSize',16)
filename = [save_dir_plot '/EXP_fig_' num2str(ifig)];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
% print('-depsc', fname);
print('-r300','-dpng', fname2);

ifig=ifig+1;
if rflags.PLOTS ==1;
    figure(ifig)
else
    figure('visible','off');
end
set(gca,'FontSize',12)
plot(1e6*beam_in(:,X_COORDINATE),1e6*beam_in(:,Y_COORDINATE),'.r')
grid on
xlabel('x [\mum]')
ylabel('y [\mum]')
filename = [save_dir_plot '/EXP_fig_' num2str(ifig)];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
% print('-depsc', fname);
print('-r300','-dpng', fname2);




e_x=(beam_in(:,X_COORDINATE));
e_xp=(beam_in(:,X_MOMENTUM)./(beam_in(:,S_MOMENTUM)));
e_y=(beam_in(:,Y_COORDINATE));
e_yp=(beam_in(:,Y_MOMENTUM)./(beam_in(:,S_MOMENTUM)));
e_s=beam_in(:,TIME_COORDINATE)';
e_ens=((beam_in(:,ENERGY_OF_PARTICLE)-mean(beam_in(:,ENERGY_OF_PARTICLE)))/mean(beam_in(:,ENERGY_OF_PARTICLE)))';
e_E=beam_in(:,ENERGY_OF_PARTICLE);

mean_e_Px=mean(beam_in(:,X_MOMENTUM));
mean_e_Py=mean(beam_in(:,Y_MOMENTUM));
std_e_Px=std(beam_in(:,X_MOMENTUM));
std_e_Py=std(beam_in(:,Y_MOMENTUM));





%e_emx_sqrt=sqrt(mean(e_x.^2)*mean(e_xp.^2)-mean(e_x.*e_xp)^2)
e_emx_sqrt=sqrt(mean((e_x-mean(e_x)).^2)*mean((e_xp-mean(e_xp)).^2)-mean((e_x-mean(e_x)).*(e_xp-mean(e_xp)))^2);
e_emy_sqrt=sqrt(mean((e_y-mean(e_y)).^2)*mean((e_yp-mean(e_yp)).^2)-mean((e_y-mean(e_y)).*(e_yp-mean(e_yp)))^2);

sigma_0x=std(e_x);
sigma_0y=std(e_y);

mean_x=mean(e_x);
mean_y=mean(e_y);


mean_XP=mean(e_xp);
mean_YP=mean(e_yp);


gamma=mean(e_E)/(0.511e6);
delta_gamma=std(e_E./(0.511e6));
norm_em_x=sqrt(gamma^2-1)*e_emx_sqrt;
norm_em_y=sqrt(gamma^2-1)*e_emy_sqrt;

energy_spread=std(e_E)/mean(e_E);
mean_e_E=mean(e_E);
std_e_E=std(e_E);
n_m_p=length(e_E);

save([save_dir_plot '/beam_mat.dat'],'e_emx_sqrt','e_emy_sqrt','sigma_0x','sigma_0y',...
    'norm_em_x','norm_em_y','energy_spread','mean_e_E','std_e_E','mean_e_Px','mean_e_Py','std_e_Px',...
    'std_e_Py','mean_x','mean_y','mean_XP','mean_YP','electron_angle');

f_name=[save_dir_plot '/beam_stat.txt'];
%  filename = sprintf('NEW_plot.txt',tilted_name);
fileID = fopen(f_name,'w');
fprintf(fileID,'Electron beam non norm Emittances:\n');
fprintf(fileID,'Emit X = %10.5e\n',e_emx_sqrt);
fprintf(fileID,'Emit Y = %10.5e\n',e_emy_sqrt);
fprintf(fileID,'Electron beam normalized Emittances:\n');
fprintf(fileID,'Emit X n = %10.5e\n',norm_em_x);
fprintf(fileID,'Emit Y n = %10.5e\n',norm_em_y);
fprintf(fileID,'\n');
fprintf(fileID,'sigma_x = %10.5e [um]\n',sigma_0x*1e6);
fprintf(fileID,'sigma_y = %10.5e [um]\n',sigma_0y*1e6);
fprintf(fileID,'mean_energy = %10.5e [MeV]\n',mean(e_E)/1e6);
fprintf(fileID,'std_energy = %10.5e [MeV]\n',std(e_E)/1e6);
fprintf(fileID,'energy_spread = %10.5e\n',energy_spread);
fprintf(fileID,'gamma = %10.5e\n',gamma);
fprintf(fileID,'delta_gamma = %10.5e\n',delta_gamma);
fprintf(fileID,'----------------------------------------------------\n');
fprintf(fileID,'NMP = %10.5e\n',n_m_p);
fprintf(fileID,'N P = %10.5e\n',n_m_p*beam_in(1,3));
fprintf(fileID,'zeta mean = %10.5e\n',mean(electron_angle));
fprintf(fileID,'zeta std = %10.5e\n',std(electron_angle));
fclose(fileID);




f_name=[save_dir_plot '/beam_stat_TEX.txt'];
%  filename = sprintf('NEW_plot.txt',tilted_name);
fileID = fopen(f_name,'w');
fprintf(fileID,'Electron beam non norm Emittances:\\\\ \n');
fprintf(fileID,'Emit X $= %10.5e $\\\\ \n',e_emx_sqrt);
fprintf(fileID,'Emit Y $= %10.5e $\\\\ \n',e_emy_sqrt);
fprintf(fileID,'Electron beam normalized Emittances:\\\\ \n');
fprintf(fileID,'Emit X n $= %10.5e $\\\\ \n',norm_em_x);
fprintf(fileID,'Emit Y n $= %10.5e $\\\\ \n',norm_em_y);
%fprintf(fileID,'\n');
fprintf(fileID,'$\\sigma_x = %10.5e [\\mu m]$\\\\ \n',sigma_0x*1e6);
fprintf(fileID,'$\\sigma_y = %10.5e [\\mu m]$ \\\\ \n',sigma_0y*1e6);
%fprintf(fileID,'\n');
fprintf(fileID,'mean energy $= %10.2f $[MeV]\\\\ \n',mean(e_E)/1e6);
fprintf(fileID,'std energy $= %10.2f $[MeV]\\\\ \n',std(e_E)/1e6);
fprintf(fileID,'energy spread $= %10.5e $\\\\ \n',energy_spread);
fprintf(fileID,'$\\gamma = %10.5e $\\\\ \n',gamma);
fprintf(fileID,'$\\delta\\gamma = %10.5e $\\\\ \n',delta_gamma);
%fprintf(fileID,'\n');
fprintf(fileID,'NMP $= %10.5e $\n',n_m_p);
fprintf(fileID,'N P $= %10.5e $\n',n_m_p*beam_in(1,3));
fprintf(fileID,'$\\zeta = %10.5e $\n',mean(electron_angle));
fclose(fileID);
%end
