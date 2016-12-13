function [theta_out]=rejection_method(E_cm)
global ev cs cross
 
flag=true;

while flag

r=min(cs(:,1))+(max(cs(:,1))-min(cs(:,1)))*rand(1);

Fmax=max(interp2(ev,cs,cross,E_cm,[min(cs(:,1)) max(cs(:,1))],'linear'));

y=Fmax*rand(1);

Vq=interp2(ev,cs,cross,E_cm,r,'linear');

if y<Vq
flag=false;     
    theta_out=r;
end

    
end