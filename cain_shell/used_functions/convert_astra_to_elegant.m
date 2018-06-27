function elegant_beam=convert_astra_to_elegant(astra_beam)
% convert astra to elegant

dim_as=size(astra_beam);
elegant_beam=zeros(dim_as(1),6);

%% convert P_z -> P_tot
reference_Pz=astra_beam(1,6);
astra_beam(1,6)=0;
elegant_beam_tempt_Pz=astra_beam(:,6)+reference_Pz;
%% create for 6th coloun is P_tot MeV
elegant_beam(:,6)=sqrt(elegant_beam_tempt_Pz.^2+astra_beam(:,4).^2+astra_beam(:,5).^2)./(0.511.*1e6);

%% remove first refefence particle
% kill Marchello reference z0 in astra is ditance from the gun for cain we take beam at zero of IP point 
reference_z=0;%astra_beam(1,3)
astra_beam(1,3)=0;
%%  create for 5th coloun is z [s]
elegant_beam(:,5)=astra_beam(:,3)./3e8+reference_z;

%% 1th coloun of elegant x coordinate
elegant_beam(:,1)=astra_beam(:,1);
%% 2th coloun of elegant xp
elegant_beam(:,2)=astra_beam(:,4)./elegant_beam_tempt_Pz;
%% 3th coloun of elegant y coordinate
elegant_beam(:,3)=astra_beam(:,2);
%% 4th coloun of elegant yp coordinate
elegant_beam(:,4)=astra_beam(:,5)./elegant_beam_tempt_Pz;


%% just plots
% for ni=1:1:6
% figure(ni)
% set(gca,'FontSize',18)
% histogram(elegant_beam(:,ni),20)
% title({['coloun ' num2str(ni) ];['std ' num2str(std(elegant_beam(:,ni))) ]...
%     ;['mean ' num2str(mean(elegant_beam(:,ni))) ]})
% 
% end




