clear all; close all; clc;

ou=['beam_B_defocusing']

out_folder_1=['_1'];
out_folder_2=['_2'];
out_folder_for_plots=['comp_B_1_2'];
mkdir([pwd '/' out_folder_for_plots '/'])

ifig=5001;
loop_number=200;

%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss

all_data_phot_1=[];
full_spectrum_1=[];
phot_angle_1=[];
weigth_1=0;
x_phot_1=[];
y_phot_1=[];
xp_phot_1=[];
yp_phot_1=[];
all_data_phot_2=[];
full_spectrum_2=[];
phot_angle_2=[];
weigth_2=0;
x_phot_2=[];
y_phot_2=[];
xp_phot_2=[];
yp_phot_2=[];

for ni=1:1:loop_number
    photons_data_1=dlmread([pwd '/' [ou out_folder_1] '/cain_tmp/cain_output_photons_' num2str(ni) '.dat'],'',1,0);%read electrons from cain
    
    x_phot_1=[x_phot_1;photons_data_1(:,5)];
    y_phot_1=[y_phot_1;photons_data_1(:,6)];
    xp_phot_1=[xp_phot_1;photons_data_1(:,9)];
    yp_phot_1=[yp_phot_1;photons_data_1(:,10)];
    weigth_1=photons_data_1(1,3);  
    
    full_spectrum_1=[full_spectrum_1;photons_data_1(:,8)./1e3];
    phot_angle_1=[phot_angle_1;atan(sqrt(photons_data_1(:,9).^2+photons_data_1(:,10).^2)./photons_data_1(:,11))];
	
	photons_data_2=dlmread([pwd '/' [ou out_folder_2] '/cain_tmp/cain_output_photons_' num2str(ni) '.dat'],'',1,0);%read electrons from cain
    
    x_phot_2=[x_phot_2;photons_data_2(:,5)];
    y_phot_2=[y_phot_2;photons_data_2(:,6)];
    xp_phot_2=[xp_phot_2;photons_data_2(:,9)];
    yp_phot_2=[yp_phot_2;photons_data_2(:,10)];    
    weigth_2=photons_data_2(1,3);   
    
    full_spectrum_2=[full_spectrum_2;photons_data_2(:,8)./1e3];
    phot_angle_2=[phot_angle_2;atan(sqrt(photons_data_2(:,9).^2+photons_data_2(:,10).^2)./photons_data_2(:,11))];
    
end

number_of_photons_1=length(full_spectrum_1)*weigth_1*weigth_1;
number_of_photons_2=length(full_spectrum_2)*weigth_2*weigth_2;
el_angel=4e-5
aa_1=find(abs(phot_angle_1)<el_angel);
aa_2=find(abs(phot_angle_2)<el_angel);



for zero_x_axes=[0 2500 4000]
nbin_plot=30;
figure(ifig)
ifig=ifig+1;
hold on
plot(linspace(zero_x_axes,max(full_spectrum_1(aa_1)),nbin_plot),smooth(hist(full_spectrum_1(aa_1),nbin_plot)*weigth_1),'-r','LineWidth',2)
plot(linspace(zero_x_axes,max(full_spectrum_2(aa_2)),nbin_plot),smooth(hist(full_spectrum_2(aa_2),nbin_plot)*weigth_2),'--b','LineWidth',2)
hold off
grid on
 ylim([zero_x_axes (max(smooth(hist(full_spectrum_1(aa_1),nbin_plot)*weigth_1))+0.15*max(smooth(hist(full_spectrum_1(aa_1),nbin_plot)*weigth_1)))])
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')

suptitle({['For theta<' num2str(el_angel) ' [rad]'];[ regexprep(out_folder_1,'_',' ') ' Total number of photons=' num2str(number_of_photons_1/loop_number,'%10.2e') ];...
[ regexprep(out_folder_1,'_',' ') ' number scattered photons in theta<' num2str(length(full_spectrum_1(aa_1))*weigth_1/loop_number,'%10.2e') ]...
;[ regexprep(out_folder_2,'_',' ') ' Total number of photons=' num2str(number_of_photons_2/loop_number,'%10.2e') ];...
[ regexprep(out_folder_2,'_',' ') ' number scattered photons in theta<' num2str(length(full_spectrum_2(aa_2))*weigth_2/loop_number,'%10.2e') ]})

% text(zero_x_axes,(max(smooth(hist(full_spectrum_2(aa_2),nbin_plot)*weigth_1))-0.1*max(smooth(hist(full_spectrum_2(aa_2),nbin_plot)*weigth_1)))...
%     ,{[ regexprep(out_folder_1,'_',' ') ' bandwidth=' num2str(std(full_spectrum_1(aa_1))/mean(full_spectrum_1(aa_1)))];...
%     [ regexprep(out_folder_2,'_',' ') ' bandwidth=' num2str(std(full_spectrum_2(aa_2))/mean(full_spectrum_2(aa_2)))]})
legend([ regexprep(out_folder_1,'_',' ') ' bandwidth=' num2str(std(full_spectrum_1(aa_1))/mean(full_spectrum_1(aa_1)))],...
    [ regexprep(out_folder_2,'_',' ') ' bandwidth=' num2str(std(full_spectrum_2(aa_2))/mean(full_spectrum_2(aa_2)))],0)

filename = [pwd '/' out_folder_for_plots '/com_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);

end


