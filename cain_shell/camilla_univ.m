
%% after propagation


bu_factor=50
theta=bu_factor*5e-5;
num_s =3;
w_size=bu_factor*5e-5;
nx_bin=50;
ny_bin=50;
xcor_allb=[]; ycor_allb=[];
    eval(['x1_coor' int2str(num_s) '=[];']);
    eval(['y1_coor' int2str(num_s) '=[];']);
eval(['aa_' num2str(num_s) '=find(abs(phot_angle_' num2str(num_s) ')<' num2str(theta) ');']);

for l=1%[0:0.0001:0.01]% 0.1:0.1:0.5];%
    close all
    
    eval(['ycor_nop=y_coor' int2str(num_s) '(aa_' int2str(num_s) ');']);
    eval(['xcor_nop=x_coor' int2str(num_s) '(aa_' int2str(num_s) ');']);
    eval(['Sz=Sz_' int2str(num_s) ';']);
    eval(['Sz_sep=Sz_' int2str(num_s) '(aa_' int2str(num_s) ');']);
    eval(['x1_coor' int2str(num_s) '=x_coor' int2str(num_s) '+(xp_' int2str(num_s) './zp_' int2str(num_s) ').*l;']);
    eval(['y1_coor' int2str(num_s) '=y_coor' int2str(num_s) '+(yp_' int2str(num_s) './zp_' int2str(num_s) ').*l;']);
    eval(['xcor_allb=x1_coor' int2str(num_s) ';']);
    eval(['ycor_allb=y1_coor' int2str(num_s) ';']);
    eval(['xcor_n=xcor_nop+(xp_' int2str(num_s) '(aa_' int2str(num_s) ')./zp_' int2str(num_s) '(aa_' int2str(num_s) ')).*l;']);
    eval(['ycor_n=ycor_nop+(yp_' int2str(num_s) '(aa_' int2str(num_s) ')./zp_' int2str(num_s) '(aa_' int2str(num_s) ')).*l;']);
    
    d1=zeros(ny_bin,nx_bin);
    d2=zeros(ny_bin,nx_bin);
    f=zeros(ny_bin,nx_bin);
    d1_theta=zeros(ny_bin,nx_bin);
    d2_theta=zeros(ny_bin,nx_bin);
    f_theta=zeros(ny_bin,nx_bin);
    
    x_l=2*w_size/nx_bin;
    y_l=2*w_size/ny_bin;
    
    
    for ni=1:1:length(xcor_allb);
        if (xcor_allb(ni)<w_size && xcor_allb(ni)>-w_size)
            if (ycor_allb(ni)<w_size && ycor_allb(ni)>-w_size)
                d1(1+floor((ycor_allb(ni)+w_size)./y_l),1+floor((xcor_allb(ni)+w_size)./x_l))=...
                    d1(1+floor((ycor_allb(ni)+w_size)./y_l),1+floor((xcor_allb(ni)+w_size)./x_l))+Sz(ni);
                d2(1+floor((ycor_allb(ni)+w_size)./y_l),1+floor((xcor_allb(ni)+w_size)./x_l))=...
                    d2(1+floor((ycor_allb(ni)+w_size)./y_l),1+floor((xcor_allb(ni)+w_size)./x_l))+1;
            end
        end
    end
    
    
    for nyi=1:1:ny_bin
        for nxi=1:1:nx_bin
            if(d2(nyi,nxi)==0)
                f(nyi,nxi)=0;
            else
                f(nyi,nxi)=d1(nyi,nxi)./d2(nyi,nxi);
            end
        end
    end
    
    
    
    for ni=1:1:length(xcor_n);
        if (xcor_n(ni)<w_size && xcor_n(ni)>-w_size)
            if (ycor_n(ni)<w_size && ycor_n(ni)>-w_size)
                d1_theta(1+floor((ycor_n(ni)+w_size)./y_l),1+floor((xcor_n(ni)+w_size)./x_l))=d1_theta(1+floor((ycor_n(ni)+w_size)./y_l),1+floor((xcor_n(ni)+w_size)./x_l))+Sz_sep(ni);
                d2_theta(1+floor((ycor_n(ni)+w_size)./y_l),1+floor((xcor_n(ni)+w_size)./x_l))=d2_theta(1+floor((ycor_n(ni)+w_size)./y_l),1+floor((xcor_n(ni)+w_size)./x_l))+1;
            end
        end
    end
    for nyi=1:1:ny_bin
        for nxi=1:1:nx_bin
            if(d2_theta(nyi,nxi)==0)
                f_theta(nyi,nxi)=0;
            else
                f_theta(nyi,nxi)=d1_theta(nyi,nxi)./d2_theta(nyi,nxi);
            end
        end
    end
    
    
    
    
    
    
    
    
    figure(100)
    % ifig=ifig+1;
    subplot 221
    set(pcolor(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        f),'EdgeColor','none')
    colorbar('SouthOutside')
    title('prop')
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    
    subplot 222
    set(pcolor(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        d2),'EdgeColor','none')
    colorbar('SouthOutside')
    title('prop')
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    suptitle(['L=' num2str(l)])
    
    
    subplot 223
    set(pcolor(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        f_theta),'EdgeColor','none')
    colorbar('SouthOutside')
    title('prop')
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    
    subplot 224
    set(pcolor(linspace(-w_size,w_size,nx_bin),...
        linspace(-w_size,w_size,ny_bin),...
        d2_theta),'EdgeColor','none')
    colorbar('SouthOutside')
    title('prop')
    set(gca,'FontSize',16)
    xlabel('X')
    ylabel('Y')
    suptitle(['L=' num2str(l) ' m; Stokes=' num2str(STOKES(num_s,:)) ])
    
    
    filename = ['EXP_plot_' num2str(1) '_' num2str(l*1e6) ];
    fname = [ filename '.png'];
    print('-dpng', fname);
    
end












