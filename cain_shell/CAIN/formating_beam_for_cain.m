function [beam_property] = formating_beam_for_cain(beam_phasespace,turn_number)
global  BASE_DIRECTORY rflags;
% global  beam_parameters
%%%%%%%%%%%%%%%%%  formating  beam data for cain %%%%%%%%%%%%%%%%%%%%%%%%%%
%  size(beam_phasespace)
%  stop

% file_electron1 = load([directory_electrons 'ELI_NP_electron_input_720MeV.dat']);
% file_electron1 = dlmread(['exinj_pls-start_0866_20pC_match_good_dd.out.sdds.asci'],'',10,0);
% beam_phasespace=dlmread(['exinj_pls-start_0866_20pC_match_good_dd.out.sdds.asci'],'',10,0);
column11 = beam_phasespace(:,1);
column12 = beam_phasespace(:,2);
column13 = beam_phasespace(:,3);
column14 = beam_phasespace(:,4);
column15 = beam_phasespace(:,5);
column16 = beam_phasespace(:,6);


% 
% NumberMP = length(beam_phasespace(:,1));
% column11 = std(beam_phasespace(:,1)).*randn(1,NumberMP)';
% column12 = std(beam_phasespace(:,2)).*randn(1,NumberMP)';
% column13 = std(beam_phasespace(:,3)).*randn(1,NumberMP)';
% column14 = std(beam_phasespace(:,4)).*randn(1,NumberMP)';
% column15 = std(beam_phasespace(:,5)).*randn(1,NumberMP)';
% % column16 = std(beam_phasespace(:,6)).*randn(1,NumberMP)';






rest_mass = 0.511;
Pz1 = rest_mass.*column16;
Px1 = column12.*Pz1;
Py1 = column14.*Pz1;
Ptot1 = sqrt(Pz1.^2 + Px1.^2 + Py1.^2);
Etot1 = sqrt(Ptot1.^2 + rest_mass.^2);
T1 = column15 - mean(column15);

% emit = 1e6*722*sqrt(mean(column11.^2) * mean(column12.^2) - mean(column11.*column12).^2);




K0=2;
genname0 = 1;
NumberMP = length(column11);
% Ne = 1.56e9; % 250 pC
%Ne = 3.125e9; % 500 pC
% rflags.chargebunch = 250e-12;%Charge per electrons bunch [c] 
Echarge = 1.60e-19;% Charge of electron [c]
Ne = rflags.chargebunch/Echarge;% Number electrons in bunch



weight0 = Ne/NumberMP;
K = ones(1,NumberMP)*K0; % 1 colomn
genname = ones(1,NumberMP)*genname0; % 2 colomn
weight = ones(1,NumberMP)*weight0; % 3 colomn


T = 3e8*T1'; % 4 colomn
X = column11'; % 5 colomn
Y = column13'; % 6 colomn
Z =zeros(1,NumberMP); % 7 colomn


energy = 1e6*Etot1'; % 8 colomn

Pxi = 1e6*Px1'; % 9 colomn
Pyi = 1e6*Py1'; % 10 colomn
Pzi = 1e6*Pz1'; % 11 colomn

Sx = zeros(1, NumberMP);% % 12 colomn
Sy = zeros(1, NumberMP);% % 13 colomn
Sz = zeros(1, NumberMP);% % 14 colomn

beam_property=[K;genname;weight;T;X;Y;Z;energy;Pxi;Pyi;Pzi;Sx;Sy;Sz;];

% beam_initial=transpose([K;genname;weight;T;X;Y;Z;energy;Pxi;Pyi;Pzi;Sx;Sy;Sz;]);
%
%
% fid = fopen(['ELI_NP_electron_input_720MeV_cain.dat'],'w'); %!!!!!! Name of the file
%
% fprintf(fid,' %i    %i       %1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e \n',transpose(beam_initial));
% fclose(fid);



if turn_number<2 && rflags.PLOTS==1
    
    ifig=0+turn_number;
    
    
    ifig=ifig+1;
    figure(ifig)
    subplot 221
    set(gca,'FontSize',12)
    hist(Px1,50)
    xlabel('P_x [MeV]')
    title(['RMS = ' num2str(std(Px1)) ' MeV'],'FontSize',9)
    subplot 222
    set(gca,'FontSize',12)
    hist(Py1,50)
    xlabel('P_y [MeV]')
    title(['RMS = ' num2str(std(Py1)) ' MeV'],'FontSize',9)
    subplot 224
    set(gca,'FontSize',12)
    hist(Ptot1,50)
    xlabel('P_{tot} [MeV]')
    title(['RMS = ' num2str(std(Ptot1)) ' MeV'],'FontSize',9)
    subplot 223
    set(gca,'FontSize',12)
    hist(Pz1,50)
    xlabel('P_z [MeV]')
    title(['RMS = ' num2str(std(Pz1)) ' MeV '],'FontSize',9)
    filename = [BASE_DIRECTORY 'initial_beam_plot/fig_' num2str(ifig)];
    fname = [ filename '.png'];
    print('-dpng', fname);
    
    ifig=ifig+1;
    figure(ifig)
    subplot 211
    set(gca,'FontSize',14)
    hist(Ptot1,50)
    xlabel('P_{tot} [MeV]')
    title(['RMS = ' num2str(std(Ptot1)) ' MeV ' ' Mean = ' num2str(mean(Ptot1))])
    subplot 212
    set(gca,'FontSize',14)
    hist(Etot1,50)
    xlabel('E_{tot} [MeV]')
    title(['RMS = ' num2str(std(Etot1)) ' MeV ' ' Mean = ' num2str(mean(Etot1))])
    filename = [BASE_DIRECTORY 'initial_beam_plot/fig_' num2str(ifig)];
    fname = [ filename '.png'];
    print('-dpng', fname);
    
    ifig=ifig+1;
    figure(ifig)
    set(gca,'FontSize',14)
    hist(1e12*T1,50)
    xlabel('T [ps]')
    title(['RMS = ' num2str(std(1e12*T1)) ' ps' ' Mean = ' num2str(mean(1e12*T1)) ' ps'])
    filename = [BASE_DIRECTORY 'initial_beam_plot/fig_' num2str(ifig)];
    fname = [ filename '.png'];
    print('-dpng', fname);
    
    
    ifig=ifig+1;
    figure(ifig)
    subplot 221
    set(gca,'FontSize',11)
    hist(1e6*column11,50)
    xlabel('x [\mum]')
    title(['RMS = ' num2str(std(1e6*column11)) ' \mum'],'FontSize',9)
    subplot 222
    set(gca,'FontSize',11)
    hist(1e6*column13,50)
    xlabel('y [\mum]')
    title(['RMS = ' num2str(std(1e6*column13)) ' \mum'],'FontSize',9)
    subplot 223
    set(gca,'FontSize',11)
    hist(1e9*column15,50)
    xlabel('T [ns]')
    title(['RMS = ' num2str(std(1e12*column15)) ' ps'],'FontSize',9)
    subplot 224
    set(gca,'FontSize',11)
    hist(rest_mass*column16,50)
    xlabel('P_z [MeV]')
    title(['RMS = ' num2str(std(rest_mass*column16)) ' MeV ' ' Mean = ' num2str(mean(rest_mass*column16))],'FontSize',9)
    filename = [BASE_DIRECTORY 'initial_beam_plot/fig_' num2str(ifig)];
    fname = [ filename '.png'];
    print('-dpng', fname);
    
    
    ifig=ifig+1;
    figure(ifig)
    subplot 221
    set(gca,'FontSize',11)
    hist(1e6*column11,50)
    xlabel('x [\mum]')
    title(['RMS = ' num2str(std(1e6*column11)) ' \mum'],'FontSize',9)
    subplot 222
    set(gca,'FontSize',11)
    hist(1e6*column13,50)
    xlabel('y [\mum]')
    title(['RMS = ' num2str(std(1e6*column13)) ' \mum'],'FontSize',9)
    subplot 223
    set(gca,'FontSize',11)
    hist(1e12*T1,50)
    xlabel('T [ps]')
    title(['RMS = ' num2str(std(1e12*T1)) ' ps'],'FontSize',9)
    subplot 224
    set(gca,'FontSize',11)
    hist(Etot1,50)
    xlabel('E_{tot} [MeV]')
    title(['RMS = ' num2str(std(Etot1)) ' MeV ' ' Mean = ' num2str(mean(Etot1))],'FontSize',9)
    filename = [BASE_DIRECTORY 'initial_beam_plot/fig_' num2str(ifig)];
    fname = [ filename '.png'];
    print('-dpng', fname);
    
    
    
    close all;
    
    
    
    
    
    
    
end




%indc_in_accept=find(beam_phasespace(16,:)<1);
%
% WEIGHT=3;
% TIME_COORDINATE=7;%now we start use s(m) for caine
% X_COORDINATE=5;
% Y_COORDINATE=6;
% ENERGY_OF_PARTICLE=8;
% X_MOMENTUM=9;
% Y_MOMENTUM=10;
% S_MOMENTUM=11;
% % POLARISATION: 12 13 14
% N_COMPTON_HIT=15;
% TURN_LAST_COMPTON_HIT=16;
%
% beam_property(:,1)=(beam_phasespace(7,:))';
%
% beam_property(:,2)=(beam_phasespace(8,:))';
%
% beam_property(:,WEIGHT)=(beam_phasespace(9,:))';
%
% beam_property(:,4)=(beam_phasespace(10,:))';
%
% beam_property(:,X_COORDINATE)=(beam_phasespace(1,:))';
%
% beam_property(:,Y_COORDINATE)=(beam_phasespace(3,:))';
%
% beam_property(:,TIME_COORDINATE)=(beam_phasespace(5,:))';
%
% beam_property(:,ENERGY_OF_PARTICLE)=beam_parameters.initial_beam_energy.*(1+(beam_phasespace(6,:))');
%
% beam_property(:,S_MOMENTUM)=sqrt((beam_property(:,ENERGY_OF_PARTICLE).^2-beam_parameters.Emass^2)./(1+((beam_phasespace(2,:))').^2+((beam_phasespace(4,:))').^2));
%
% beam_property(:,X_MOMENTUM)=((beam_phasespace(2,:))').*(beam_property(:,S_MOMENTUM));
%
% beam_property(:,Y_MOMENTUM)=((beam_phasespace(4,:))').*(beam_property(:,S_MOMENTUM));
%
% beam_property(:,12)=(beam_phasespace(11,:))';
%
% beam_property(:,13)=(beam_phasespace(12,:))';
%
% beam_property(:,14)=(beam_phasespace(13,:))';
%
