
theta_angle=10e-6;
diapason=(10:1:101);
el_angel=0;
qq=0; num_phot_th=[]; bandwith_cm=[]; bandwith_non_norm=[];
for ni=diapason;
    qq=qq+1;
    el_angel_1=ni*theta_angle;
    aa_1=find(abs(phot_angle)<el_angel_1);
    
    num_phot_th(qq)=length(full_spectrum(aa_1))*weigth;
    bandwith_cm(qq)=std(full_spectrum(aa_1))/mean(full_spectrum(aa_1));
    
    mean_Sx(qq)=mean(Sx(aa_1));
    mean_Sy(qq)=mean(Sy(aa_1));
    mean_Sz(qq)=mean(Sz(aa_1));
    
    if bandwith_cm(qq)>=0.0048 && bandwith_cm(qq)<=0.0051 % find angle for given bandwith
        el_angel=el_angel_1
        % stop;
    end
end
aa=find(abs(phot_angle)<el_angel);

bandwith=std(full_spectrum(aa))/mean(full_spectrum(aa))
el_angel
l=10;
    xcor_n=x_phot(aa)+(xp_phot(aa)./zp_phot(aa)).*l;
    ycor_n=y_phot(aa)+(yp_phot(aa)./zp_phot(aa)).*l;

for nbin_plot=[10 20 30 50 100 1000 1e4]
xs=linspace(min(xcor_n),max(xcor_n),nbin_plot);
ys=smooth(hist(xcor_n,nbin_plot));
iysm=find(ys==max(ys));

C1 = abs(bsxfun(@minus,ys(1:iysm)',max(ys)/2))
[~,idx1] = min(C1(:,1:size(C1,2)))
C2 = abs(bsxfun(@minus,ys(iysm:length(ys))',ys(idx1)))
[~,idx2] = min(C2(:,1:size(C2,2)))
idx2=idx2+iysm-1;

fwhm_x=abs(xs(idx2)-xs(idx1))
std(xcor_n)

figure(ifig)
ifig=ifig+1;
set(gca,'FontSize',16)
hold on
hist(xcor_n,nbin_plot)
plot(xs,ys,'-r','LineWidth',3)
plot(xs(idx1),ys(idx1),'og','LineWidth',5)
% plot(xs(iysm:length(ys)),ys(iysm:length(ys)),'og','LineWidth',5)
plot(xs(idx2),ys(idx2),'sg','LineWidth',5)
hold off
grid on
xlabel('X [m]')
title({['L=' num2str(l) ' m; Stokes=' num2str(rflags.STOKES) '; nbin=' num2str(nbin_plot) ];[' BW=' num2str(bandwith,'%10.2e') ' for \theta=' num2str(el_angel,'%10.2e') ];...
    ['fwhm =' num2str(fwhm_x*1e3) ' mm; std(x)=' num2str(std(xcor_n)*1e3) ' mm;']})
filename = [ out_folder  'X_photons_plot_' num2str(ifig) ];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
%print('-depsc', fname);
print('-r300','-dpng', fname2);

end


for nbin_plot=[10 20 30 50 100 1000 1e4]
xs=linspace(min(ycor_n),max(ycor_n),nbin_plot);
ys=smooth(hist(ycor_n,nbin_plot));
iysm=find(ys==max(ys));

C1 = abs(bsxfun(@minus,ys(1:iysm)',max(ys)/2))
[~,idx1] = min(C1(:,1:size(C1,2)))
C2 = abs(bsxfun(@minus,ys(iysm:length(ys))',ys(idx1)))
[~,idx2] = min(C2(:,1:size(C2,2)))
idx2=idx2+iysm-1;

fwhm_y=abs(xs(idx2)-xs(idx1))
std(ycor_n)

figure(ifig)
ifig=ifig+1;
set(gca,'FontSize',16)
hold on
hist(ycor_n,nbin_plot)
plot(xs,ys,'-r','LineWidth',3)
plot(xs(idx1),ys(idx1),'og','LineWidth',5)
% plot(xs(iysm:length(ys)),ys(iysm:length(ys)),'og','LineWidth',5)
plot(xs(idx2),ys(idx2),'sg','LineWidth',5)
hold off
grid on
xlabel('Y [m]')
title({['L=' num2str(l) ' m; Stokes=' num2str(rflags.STOKES) '; nbin=' num2str(nbin_plot) ];[' BW=' num2str(bandwith,'%10.2e') ' for \theta=' num2str(el_angel,'%10.2e') ];...
    ['fwhm =' num2str(fwhm_y*1e3) ' mm; std(y)=' num2str(std(xcor_n)*1e3) ' mm;']})
filename = [ out_folder  'Y_photons_plot_' num2str(ifig) ];
fname = [ filename '.eps'];fname2 = [ filename '.png'];
% print('-depsc', fname);
print('-r300','-dpng', fname2);

end




