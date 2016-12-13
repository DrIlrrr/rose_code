function load_crossection_for_gamma_gamma
% function [ev,cs,cross]=load_crossection_for_gamma_gamma
global rflags ev cs cross

if rflags.gamma_gamma==1;
    %     ev=load('used_functions/e_grid.dat','ascii')';
    %     cs=load('used_functions/w_grid.dat','ascii')';
    %     cross=load('used_functions/cs_grid.dat','ascii')';
    
    ev=load('used_functions/rs_grid.dat','ascii')';
    cs=load('used_functions/w_grid.dat','ascii')';
    cross=load('used_functions/cs_grid.dat','ascii')';
    
%     cs_temp=acos(cs(:,1));
%     cs=[];
%     cs=cs_temp;
    
    plot_cross=0;
    if plot_cross==1;
        ev_m=linspace(0.1,2.4,1001);
        cthv_m=linspace(-0.9995,0.9995,101);
                
        [Eq,Wq] = meshgrid(ev_m,cthv_m');
                
        % Eq=ev;
        % Wq=cs;
        
        Vq=interp2(ev,cs,cross,Eq,Wq,'linear');
        
        
        
        figure(90)
        subplot(2,2,1)
        surf(Wq,Eq,Vq)
        title('interpol')
        subplot(2,2,2)
%         mesh(acos(cs(:,1)),ev(1,:),cross')
        mesh((cs(:,1)),ev(1,:),cross')
        title('data')
        subplot(2,2,[3 4])
        hold on
        surf(Wq,Eq,Vq)
        %  scatter3(Wq,Eq,Vq,'r')
        % scatter3(grid300(:,1),grid300(:,2),grid300(:,3),'.b')
        mesh(cs(:,1),ev(1,:),cross')%,'EdgeColor','none')
        hold off
%         xlabel('Cos(\theta)')
        xlabel('\theta')
        ylabel('E_{cm}')
        view(44,30)
        
        
    end
    
end