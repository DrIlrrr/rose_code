function [theta_out,cross]=rejection_method_approx_gamma_gamma(E_cm)
% global ev cs cross
stop
flag=true;

while flag
    
    
     r=-1+2*rand(1);
    
    %     Fmax=max(interp2(ev,cs,cross,E_cm,[min(cs(:,1)) max(cs(:,1))],'linear'));
    Fmax=max(aprox_gamma_gamma_crossection_in_the_CoM(E_cm,-1:1e-3:1));
    y=Fmax*rand(1);
    
    % Vq is value of differential crossection taked from interpolated funtion
    %     Vq=interp2(ev,cs,cross,E_cm,r,'linear');
    Vq=aprox_gamma_gamma_crossection_in_the_CoM(E_cm,r);
    if y<Vq
        flag=false;
        theta_out=r;
    end
    
%disp(['flag = ' num2str(flag)])

end


cross=Vq.*1e34;