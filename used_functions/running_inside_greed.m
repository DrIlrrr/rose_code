function [p_cm_step,gamma_cm_step,Ecm_step,NumP_1,NumP_2,delta_x,delta_y,delta_z,pair_step]=running_inside_greed(beam_1,beam_2,x_max,x_min,nx_bin,y_max,y_min,ny_bin,z_max,z_min,nz_bin)
global rflags ifig

% global ev cs cross
% we can put here findind_min_max_of_grid(beam_1,beam_2); it will make a
% briathing cell
NumP_1=zeros(nx_bin,ny_bin,nz_bin);
NumP_2=zeros(nx_bin,ny_bin,nz_bin);
Ecm_step=[];gamma_cm_step=[];p_cm_step=[];pair_step=[];

x_1=beam_1(:,5);
y_1=beam_1(:,6);
z_1=beam_1(:,4);
E_1=beam_1(:,8);


x_2=beam_2(:,5);
y_2=beam_2(:,6);
z_2=beam_2(:,4);
E_2=beam_2(:,8);

delta_x=(x_max-x_min)/nx_bin;
delta_y=(y_max-y_min)/ny_bin;
delta_z=(z_max-z_min)/nz_bin;

% Ecm_incide=[];gamma_cm_incide_cell=[];


for niy=1:1:ny_bin%+1
    for nix=1:1:nx_bin%+1
        for niz=1:1:nz_bin%+1
            
            za_1=[];za_2=[];xa_1=[];xa_2=[];ya_1=[];ya_2=[];
            
            
            za_1=find((z_1>=z_min+(niz-1)*delta_z & z_1<=z_min+(niz)*delta_z));
            za_2=find((z_2>=z_min+(niz-1)*delta_z & z_2<=z_min+(niz)*delta_z));
            
            xa_1=find((x_1(za_1)>=x_min+(nix-1)*delta_x & x_1(za_1)<=x_min+(nix)*delta_x));
            xa_2=find((x_2(za_2)>=x_min+(nix-1)*delta_x & x_2(za_2)<=x_min+(nix)*delta_x));
            
            
            ya_1=find((y_1(xa_1)>=y_min+(niy-1)*delta_y & y_1(xa_1)<=y_min+(niy)*delta_y));
            ya_2=find((y_2(xa_2)>=y_min+(niy-1)*delta_y & y_2(xa_2)<=y_min+(niy)*delta_y));
            
            beam_1_in_cell=[];
            beam_2_in_cell=[];
            
            if length(za_1)>0 &&  length(za_2)>0 && length(xa_1)>0 && length(xa_2)>0 && length(ya_1)>0 &&  length(ya_2)>0
                % this if check that in the one cell will particle of the
                % both type
                
                beam_1_in_cell=beam_1(ya_1,:);
                beam_2_in_cell=beam_2(ya_2,:);
                
                %% this is just for debuging. is show cell coordinate and number of particle both type on it
                %                 if rflags.PLOTS ==1;
                %                     fprintf('cell %4.0f %4.0f %4.0f containing P1 %4.0f P2 %4.0f \n',nix,niy,niz,length(ya_1),length(ya_2))
                %                 end
                Ecm_pair=[];gamma_cm_pair=[];p_cm_pair=[];qq=0;crosssection_pair=[];
                %% sort out of all posible pair of particle of two type inside cell
                for goi1=1:1:length(ya_1)
                    for goi2=1:1:length(ya_2)
                        qq=qq+1;
                        work_BW=0;
                        %  1  2         3     4    5    6    7     8      9        10       11    12 13 14
                        %  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss
                        
                        %% monte carlo event generator
                        cos_alpha=[];cos_theta=[];theta_3=[];theta_4=[];phi_3=[];phi_4=[];E_3=[];E_4=[];Vq=[];Ecm_pair=[];gamma_cm_pair=[];
                        %% Gamma-Gamma
                        if rflags.gamma_gamma==1;
%                             m_1=0;% photon rest mass
%                             m_2=0;% photon rest mass
%                             Ecm_pair=sqrt(m_1^2+m_2^2+2*beam_1_in_cell(goi1,8)*beam_2_in_cell(goi2,8)-2*(beam_1_in_cell(goi1,9)*beam_2_in_cell(goi2,9)+beam_1_in_cell(goi1,10)*beam_2_in_cell(goi2,10)+beam_1_in_cell(goi1,11)*beam_2_in_cell(goi2,11)));
% 
%                              if Ecm_pair>0.001e6 %if smaller cross-section is too small
%                                 work_BW=0;
                                [cos_alpha,cos_theta,theta_3,theta_4,phi_3,phi_4,E_3,E_4,Vq,Ecm_pair,gamma_cm_pair]=gamma_gamma_monte_carlo_event_generator...
                                    (beam_1_in_cell(goi1,8),beam_1_in_cell(goi1,9),beam_1_in_cell(goi1,10),beam_1_in_cell(goi1,11),beam_1_in_cell(goi1,12),beam_1_in_cell(goi1,13),beam_1_in_cell(goi1,14),...
                                    beam_2_in_cell(goi2,8),beam_2_in_cell(goi2,9),beam_2_in_cell(goi2,10),beam_2_in_cell(goi2,11),beam_2_in_cell(goi2,12),beam_2_in_cell(goi2,13),beam_2_in_cell(goi2,14));
%                                 else
%                                 work_BW=1;
%                             end
                            %% Compton Back Scattering
                        elseif rflags.compton==1;
                            [cos_alpha,cos_theta,theta_3,theta_4,phi_3,phi_4,E_3,E_4,Vq,Ecm_pair,gamma_cm_pair]=compton_monte_carlo_event_generator...
                                (beam_1_in_cell(goi1,8),beam_1_in_cell(goi1,9),beam_1_in_cell(goi1,10),beam_1_in_cell(goi1,11),beam_1_in_cell(goi1,12),beam_1_in_cell(goi1,13),beam_1_in_cell(goi1,14),...
                                beam_2_in_cell(goi2,8),beam_2_in_cell(goi2,9),beam_2_in_cell(goi2,10),beam_2_in_cell(goi2,11),beam_2_in_cell(goi2,12),beam_2_in_cell(goi2,13),beam_2_in_cell(goi2,14));
                            %% BW breit_wheeler
                        elseif rflags.breit_wheeler==1;
                            m_1=0;
                            m_2=0;%photon rest mass
                            Ecm_pair=sqrt(m_1^2+m_2^2+2*beam_1_in_cell(goi1,8)*beam_2_in_cell(goi2,8)-2*(beam_1_in_cell(goi1,9)*beam_2_in_cell(goi2,9)+beam_1_in_cell(goi1,10)*beam_2_in_cell(goi2,10)+beam_1_in_cell(goi1,11)*beam_2_in_cell(goi2,11)));
                            
                            if Ecm_pair>1.022e6 %if smaller cross-section is 0
                                work_BW=0;
                                [cos_alpha,cos_theta,theta_3,theta_4,phi_3,phi_4,E_3,E_4,Vq,Ecm_pair,gamma_cm_pair]=breit_wheeler_monte_carlo_event_generator...
                                    (beam_1_in_cell(goi1,8),beam_1_in_cell(goi1,9),beam_1_in_cell(goi1,10),beam_1_in_cell(goi1,11),beam_1_in_cell(goi1,12),beam_1_in_cell(goi1,13),beam_1_in_cell(goi1,14),...
                                    beam_2_in_cell(goi2,8),beam_2_in_cell(goi2,9),beam_2_in_cell(goi2,10),beam_2_in_cell(goi2,11),beam_2_in_cell(goi2,12),beam_2_in_cell(goi2,13),beam_2_in_cell(goi2,14));
                            else
                                work_BW=1;
                            end
                            %%    Muon production
                        elseif rflags.e_gamma==1;
                            m_1=0.511e6;% electron rest mass [eV]
                            m_2=0;% photon rest mass
                            Ecm_pair=sqrt(m_1^2+m_2^2+2*beam_1_in_cell(goi1,8)*beam_2_in_cell(goi2,8)-2*(beam_1_in_cell(goi1,9)*beam_2_in_cell(goi2,9)+beam_1_in_cell(goi1,10)*beam_2_in_cell(goi2,10)+beam_1_in_cell(goi1,11)*beam_2_in_cell(goi2,11)));
                            
                            if Ecm_pair>211.71e6;
                                gamma_cm_pair = (beam_1(goi1,8)+beam_2(goi2,8))./Ecm_pair;
                                
                                Ecm_step=[Ecm_step; Ecm_pair];
                                gamma_cm_step=[gamma_cm_step; gamma_cm_pair];
                                
                                Vq=0; E_3=0; E_4=0; theta_3=0; theta_4=0; phi_3=0; phi_4=0; cos_alpha=0; cos_theta=0;
                            end
                            %%   TPP
                        elseif rflags.TPP==1;
                            m_1=0.511e6;% electron rest mass [eV]
                            m_2=0;% photon rest mass
                            Ecm_pair=sqrt(m_1^2+m_2^2+2*beam_1_in_cell(goi1,8)*beam_2_in_cell(goi2,8)-2*(beam_1_in_cell(goi1,9)*beam_2_in_cell(goi2,9)+beam_1_in_cell(goi1,10)*beam_2_in_cell(goi2,10)+beam_1_in_cell(goi1,11)*beam_2_in_cell(goi2,11)));
                            
                            if Ecm_pair>0;
                                gamma_cm_pair = (beam_1(goi1,8)+beam_2(goi2,8))./Ecm_pair;
                                
                                Ecm_step=[Ecm_step; Ecm_pair];
                                gamma_cm_step=[gamma_cm_step; gamma_cm_pair];
                                
                                Vq=0; E_3=0; E_4=0; theta_3=0; theta_4=0; phi_3=0; phi_4=0; cos_alpha=0; cos_theta=0;
                            end
                            %% Moller
                        elseif rflags.moller==1;
                            m_1=0.511e6;% electron rest mass [eV]
                            m_2=0.511e6;% electron rest mass [eV]
                            
                            %Ecm_pair=sqrt(m_1^2+m_2^2+2*beam_1_in_cell(goi1,8)*beam_2_in_cell(goi2,8)-2*(beam_1_in_cell(goi1,9)*beam_2_in_cell(goi2,9)+beam_1_in_cell(goi1,10)*beam_2_in_cell(goi2,10)+beam_1_in_cell(goi1,11)*beam_2_in_cell(goi2,11)));
                            %gamma_cm_pair = (beam_1(goi1,8)+beam_2(goi2,8))./Ecm_pair;
                            
                            [cos_alpha,cos_theta,theta_3,theta_4,phi_3,phi_4,E_3,E_4,Vq,Ecm_pair,gamma_cm_pair]=moller_monte_carlo_event_generator...
                                (beam_1_in_cell(goi1,8),beam_1_in_cell(goi1,9),beam_1_in_cell(goi1,10),beam_1_in_cell(goi1,11),beam_1_in_cell(goi1,12),beam_1_in_cell(goi1,13),beam_1_in_cell(goi1,14),...
                                beam_2_in_cell(goi2,8),beam_2_in_cell(goi2,9),beam_2_in_cell(goi2,10),beam_2_in_cell(goi2,11),beam_2_in_cell(goi2,12),beam_2_in_cell(goi2,13),beam_2_in_cell(goi2,14));
                            %  may be this must be uncoment due to the  work_BW=0;
                            %  try if Moller will crush
                            %                                Ecm_step=[Ecm_step; Ecm_pair];
                            %                                gamma_cm_step=[gamma_cm_step; gamma_cm_pair];
                            
                        end
                        
                        if (rflags.e_gamma+rflags.TPP)==0;
                            
                            if  work_BW==0;
                                
                                pair_info=[Vq E_3 E_4 theta_3 theta_4 phi_3 phi_4 cos_alpha cos_theta Ecm_pair gamma_cm_pair];
                                pair_step=[pair_info; pair_step];
                                
                                Ecm_step=[Ecm_step; Ecm_pair];
                                gamma_cm_step=[gamma_cm_step; gamma_cm_pair];
                            end
                        end
                        
                        
                        
                    end
                end
                
                
                
                %                 NumP_1(nix,niy,niz)=length(ya_1)*beam_1(2,3);
                %                 NumP_2(nix,niy,niz)=length(ya_2)*beam_2(2,3);
                %                 local_luminosity(nix,niy,niz)=av_cross.*NumP_1(nix,niy,niz).*NumP_2(nix,niy,niz)./(delta_x.*delta_y);
                
                
            else
                %                 stop
                %                 NumP_1(nix,niy,niz)=0;
                %                 NumP_2(nix,niy,niz)=0;
                %                 local_luminosity(nix,niy,niz)=0;
                
            end
            
            
            
        end
    end
end
