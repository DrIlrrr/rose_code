%%%%%%%%%%%%%%%%DO NOT CHANGE ANYTHING BELOW (AFTER THIS LINE)%%%%%%%%%%%
% Creating the electron beam

function[beam_phasespace] = electron_beam_initial
global BASE_DIRECTORY beam_parameters home_dir;
for_plot=0;
ifig=1;
randn_for_all=randn(1,beam_parameters.NUMBER_OF_MACROPARTICLES);
randn_for_all=randn_for_all-mean(randn_for_all);
bunch_length_NO_OFFSET = beam_parameters.SPEED_OF_LIGHT*beam_parameters.bunch_length_initial*randn(1,beam_parameters.NUMBER_OF_MACROPARTICLES); % vector of initial bunch length [m]
bunch_length_NO_OFFSET = bunch_length_NO_OFFSET -mean(bunch_length_NO_OFFSET);

if(beam_parameters.BANANA==1)
if(beam_parameters.s_positions_offset==0)
    if(beam_parameters.phase_offset~=0)
        beam_parameters.s_positions_offset=beam_parameters.phase_offset/((2*pi*beam_parameters.harmonic_number/(beam_parameters.period_revol))*(180/pi)/beam_parameters.SPEED_OF_LIGHT); % change mean S position [m]
    end
end

if(for_plot==1);
    figure(ifig)
    ifig=ifig+1;
    hist(bunch_length_NO_OFFSET)
    title({['s mean no offset' num2str(mean(bunch_length_NO_OFFSET)) ] ; ['s std no offset' num2str(std(bunch_length_NO_OFFSET)) ] ; ['t mean no offset' num2str(mean(bunch_length_NO_OFFSET/beam_parameters.SPEED_OF_LIGHT)) ]})
    filename = ['plot_initial_dist_' num2str(ifig)];
    fname = [BASE_DIRECTORY filename '.png'];
    print('-dpng', fname);
end

%     energy_deviation = beam_parameters.initial_beam_energy+beam_parameters.initial_beam_energy*beam_parameters.energy_spread_initial*randn_for_all;% vector of initial electron beam energy

%     energy_deviation = beam_parameters.initial_beam_energy.*(1+beam_parameters.energy_spread_initial.*randn_for_all).*cos((3.*1e9.*bunch_length)./beam_parameters.SPEED_OF_LIGHT);% vector of initial electron beam energy


LINAC_frequency=3*1e9;
compencation_for_banana=beam_parameters.initial_beam_energy-mean(beam_parameters.initial_beam_energy.*(1+beam_parameters.energy_spread_initial.*randn_for_all).*cos((2.*pi.*LINAC_frequency.*bunch_length_NO_OFFSET)./beam_parameters.SPEED_OF_LIGHT));
energy_deviation_NO_OFFSET = beam_parameters.initial_beam_energy.*(1+beam_parameters.energy_spread_initial.*randn_for_all).*cos((2.*pi.*LINAC_frequency.*bunch_length_NO_OFFSET)./beam_parameters.SPEED_OF_LIGHT)+compencation_for_banana*beam_parameters.COMP;% vector of initial electron beam energy
energy_deviation = energy_deviation_NO_OFFSET+beam_parameters.energy_offset;  % add linac offse
bunch_length =bunch_length_NO_OFFSET+beam_parameters.s_positions_offset;  % add linac offset
%energy_deviation = beam_parameters.initial_beam_energy.*(1+beam_parameters.energy_spread_initial.*randn_for_all);% vector of initial electron beam energy
if(for_plot==1);
    
    figure(ifig)
    ifig=ifig+1;
    hist(bunch_length)
    title({['s mean' num2str(mean(bunch_length)) ] ; ['s std' num2str(std(bunch_length)) ] ; ['t mean' num2str(mean(bunch_length/beam_parameters.SPEED_OF_LIGHT)) ]})
    filename = ['plot_initial_dist_' num2str(ifig)];
    fname = [BASE_DIRECTORY filename '.png'];
    print('-dpng', fname);
    
    
    figure(ifig)
    ifig=ifig+1;
    hist(energy_deviation)
    title(['energy mean' num2str(mean(energy_deviation))])
    filename = ['initial_mean_energy_spread'];
    fname = [BASE_DIRECTORY filename '.png'];
    print('-dpng', fname);
    filename = ['plot_initial_' num2str(ifig)];
    fname = [BASE_DIRECTORY filename '.png'];
    print('-dpng', fname);
    
    figure(ifig)
    ifig=ifig+1;
    plot((bunch_length),energy_deviation,'.r')
    grid on
    title({['energy offset ' num2str(beam_parameters.energy_offset) '[eV]' ] ;
        ['S offset ' num2str(beam_parameters.s_positions_offset) '[m]' ]  ;
        ['std s ' num2str(std(bunch_length)) ];
        ['s mean ' num2str(mean(bunch_length)) ];
        ['s mean before offset ' num2str(mean(bunch_length_NO_OFFSET)) ];
        ['std energy ' num2str(std(energy_deviation)) ];
        ['mean energy before offset ' num2str(mean(energy_deviation_NO_OFFSET)) ];
        ['mean energy ' num2str(mean(energy_deviation)) ]})
    xlabel('s [m]')
    ylabel('energy [eV]')
    filename = ['plot_initial_' num2str(ifig)];
    fname = [BASE_DIRECTORY filename '.png'];
    print('-dpng', fname);
    
    
    figure(ifig)
    ifig=ifig+1;
    plot((bunch_length_NO_OFFSET),energy_deviation_NO_OFFSET,'.g')
    hold on
    plot((bunch_length),energy_deviation,'.r')
    hold off
    grid on
    set(gca,'FontSize',14)
    title({['energy offset ' num2str(beam_parameters.energy_offset) '[eV]' ] ;
        ['S offset ' num2str(beam_parameters.s_positions_offset) '[m]' ]  ;
        ['std s ' num2str(std(bunch_length)) ];
        ['s mean ' num2str(mean(bunch_length)) ];
        ['s mean before offset ' num2str(mean(bunch_length_NO_OFFSET)) ];
        ['std energy ' num2str(std(energy_deviation)) ];
        ['mean energy before offset ' num2str(mean(energy_deviation_NO_OFFSET)) ];
        ['mean energy ' num2str(mean(energy_deviation)) ]})
    legend('before offset','with offset')
    set(gca,'FontSize',16)
    xlabel('s [m]')
    ylabel('energy [eV]')
    filename = ['plot_initial_' num2str(ifig)];
    fname = [BASE_DIRECTORY filename '.png'];
    print('-dpng', fname);
    
    
    ifig=ifig+1;
    xedges = linspace(-5*1e-3,5*1e-3,1e2); yedges = linspace(5.35*1e7,5.55*1e7,1e2);
    histmat = hist2(bunch_length, energy_deviation, xedges, yedges);
    figure;
    
    set(pcolor(xedges,yedges,histmat'),'EdgeColor','none')
    shading interp
    colorbar
    set(gca,'FontSize',14)
    title({['energy offset ' num2str(beam_parameters.energy_offset) '[eV]' ] ;
        ['S offset ' num2str(beam_parameters.s_positions_offset) '[m]' ]  ;
        ['std s ' num2str(std(bunch_length)) ];
        ['s mean ' num2str(mean(bunch_length)) ];
        ['s mean before offset ' num2str(mean(bunch_length_NO_OFFSET)) ];
        ['std energy ' num2str(std(energy_deviation)) ];
        ['mean energy before offset ' num2str(mean(energy_deviation_NO_OFFSET)) ];
        ['mean energy ' num2str(mean(energy_deviation)) ]})
    set(gca,'FontSize',16)
    xlabel('S [m]');
    ylabel('relative energy [eV]');
    filename = ['plot_initial_' num2str(ifig)];
    fname = [BASE_DIRECTORY filename '.png'];
    print('-dpng', fname);
    
    
end
else
    
energy_deviation= beam_parameters.initial_beam_energy.*(1+beam_parameters.energy_spread_initial.*randn_for_all);% vector of initial electron beam energy
energy_deviation = energy_deviation+beam_parameters.energy_offset;  % add linac offse   

bunch_length_NO_OFFSET = beam_parameters.SPEED_OF_LIGHT*beam_parameters.bunch_length_initial*randn(1,beam_parameters.NUMBER_OF_MACROPARTICLES); % vector of initial bunch length [m]
bunch_length_NO_OFFSET = bunch_length_NO_OFFSET -mean(bunch_length_NO_OFFSET);

if(beam_parameters.s_positions_offset==0)
    if(beam_parameters.phase_offset~=0)
        beam_parameters.s_positions_offset=beam_parameters.phase_offset/((2*pi*beam_parameters.harmonic_number/(beam_parameters.period_revol))*(180/pi)/beam_parameters.SPEED_OF_LIGHT); % change mean S position [m]
    end
end
bunch_length =bunch_length_NO_OFFSET+beam_parameters.s_positions_offset;  % add linac offset
end

%%%%%%%%%%%%%%%%%%%%%%%


Xdeviation =beam_parameters.sigma_e_x*randn(1,beam_parameters.NUMBER_OF_MACROPARTICLES)+beam_parameters.x_positions_offset;% vector of initial vertical electron beam size m
Ydeviation =beam_parameters.sigma_e_y*randn(1,beam_parameters.NUMBER_OF_MACROPARTICLES)+beam_parameters.x_positions_offset;% vector of initial horizontal electron beam size m


% std(Xdeviation)
% stop
%data for 2 first columns in cain standart data file
K=2;% for cain data output
genname=1;% for cain data output
K0 = ones(1,beam_parameters.NUMBER_OF_MACROPARTICLES)*K;
genname0 = ones(1,beam_parameters.NUMBER_OF_MACROPARTICLES)*genname;
weight0 = ones(1,beam_parameters.NUMBER_OF_MACROPARTICLES)*beam_parameters.weight;
t =zeros(1,beam_parameters.NUMBER_OF_MACROPARTICLES);% for time of work cain
% Creating momentum
momentum = sqrt(energy_deviation.^2-beam_parameters.Emass^2);%/SPEED_OF_LIGHT;normalezd for cain % Momentum full aka Pi [Mev/c]

%%%%% NOW WE DON'T NEED CREATE MOMENTUM. WE DON'T USE IT IN PHASESPACE.
%%%%% WE USE X & Y prime and we generate it using sigma and BETAX & BETAY
% momangle = rand(1,beam_parameters.NUMBER_OF_MACROPARTICLES);
X_prime=(beam_parameters.sigma_e_x/beam_parameters.Betax)*randn(1,beam_parameters.NUMBER_OF_MACROPARTICLES);%*cos(momangle);
Y_prime=(beam_parameters.sigma_e_y/beam_parameters.Betay)*randn(1,beam_parameters.NUMBER_OF_MACROPARTICLES);%*sin(momangle);

Pzi = momentum./sqrt(1+X_prime.^2+Y_prime.^2);% Longitudinal momentum[Mev/c]
Pxi = X_prime.*Pzi; % Momentum projections[Mev/c]
Pyi = Y_prime.*Pzi;% Momentum projections[Mev/c]

%a = momentum.^2-Pzi.^2-Pti.^2;% Test for momentum
% Creating polarizations
Sx = zeros(1,beam_parameters.NUMBER_OF_MACROPARTICLES);% X polarizations
Sy = zeros(1,beam_parameters.NUMBER_OF_MACROPARTICLES);% Y polarizations
Ss = zeros(1,beam_parameters.NUMBER_OF_MACROPARTICLES);% s polarizations
compton_hiting= zeros(1,beam_parameters.NUMBER_OF_MACROPARTICLES);% information for number of compton colligionscompton_hiting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beam_property_initial=transpose([K0;genname0;weight0;0*bunch_length;Xdeviation;Ydeviation;bunch_length;energy_deviation;Pxi;Pyi;Pzi;Sx;Sy;Ss;compton_hiting;compton_hiting;]);

WEIGHT=3;
TIME_COORDINATE=7;%now we start use s(m) for caine
X_COORDINATE=5;
Y_COORDINATE=6;
ENERGY_OF_PARTICLE=8;
X_MOMENTUM=9;
Y_MOMENTUM=10;
S_MOMENTUM=11;
% POLARISATION: 12 13 14
N_COMPTON_HIT=15;
TURN_LAST_COMPTON_HIT=16;


phasespace(1,:)=(beam_property_initial(:,X_COORDINATE))';

%!!!!!!!!!!!!!!!!!!!!!!!!!!!! It's give X_prime & Y_prime ONLY IN
%IP!!!!!!!!!!
%  phasespace(2,:)=(beam_property_initial(:,X_MOMENTUM)./(beam_property_initial(:,S_MOMENTUM)))';
phasespace(2,:)=X_prime;%((beam_parameters.sigma_e_x/beam_parameters.Betax)*randn(1,beam_parameters.NUMBER_OF_MACROPARTICLES))';

phasespace(3,:)=(beam_property_initial(:,Y_COORDINATE))';

%  phasespace(4,:)=(beam_property_initial(:,Y_MOMENTUM)./(beam_property_initial(:,S_MOMENTUM)))';
phasespace(4,:)=Y_prime;%((beam_parameters.sigma_e_y/beam_parameters.Betay)*randn(1,beam_parameters.NUMBER_OF_MACROPARTICLES))';

phasespace(5,:)=beam_property_initial(:,TIME_COORDINATE)';

phasespace(6,:)=((beam_property_initial(:,ENERGY_OF_PARTICLE)-beam_parameters.initial_beam_energy)/beam_parameters.initial_beam_energy)';

beam_phasespace=zeros(16,beam_parameters.NUMBER_OF_MACROPARTICLES);
beam_phasespace(1:6,:)=phasespace(1:6,:);
beam_phasespace(7,:)=K0;
beam_phasespace(8,:)=genname0;
beam_phasespace(9,:)=weight0;
beam_phasespace(10,:)=0*bunch_length;
beam_phasespace(11,:)=Sx;
beam_phasespace(12,:)=Sy;
beam_phasespace(13,:)=Ss;
beam_phasespace(14,:)=compton_hiting;
beam_phasespace(15,:)=compton_hiting;
beam_phasespace(16,:)=zeros(1,beam_parameters.NUMBER_OF_MACROPARTICLES);% info for acceptance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([BASE_DIRECTORY '/initial_beam_data.dat'],'beam_phasespace');

if(for_plot==1);
    figure(ifig)
    ifig=ifig+1;
    hist(phasespace(5,:))
    title(['s mean ' num2str(mean(phasespace(5,:))) ])
    filename = ['plot_initial_' num2str(ifig)];
    fname = [BASE_DIRECTORY filename '.png'];
    print('-dpng', fname);
    
    
    figure(ifig)
    ifig=ifig+1;
    plot((phasespace(5,:)),(phasespace(6,:)),'.r')
    grid on
    title({['s mean ' num2str(mean(phasespace(5,:))) ];['std relative energy ' num2str(std(phasespace(6,:))) ];['mean relative energy ' num2str(mean(phasespace(6,:))) ]})
    xlabel('s [m]')
    ylabel('relative energy (E_{initial}-E_i)/E_{initial}')
    filename = ['initial_mean_S'];
    fname = [BASE_DIRECTORY filename '.png'];
    print('-dpng', fname);
    filename = ['plot_initial_' num2str(ifig)];
    fname = [BASE_DIRECTORY filename '.png'];
    print('-dpng', fname);
     stop
end


%
%
% std_S=std((phasespace(5,:)));
% std_relative_en=std((phasespace(6,:)));
%
% std_X=std(phasespace(1,:));
% std_X_prime=std(phasespace(2,:));
% std_Y_prime=std(phasespace(4,:));