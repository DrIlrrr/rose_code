function plots_beam(beam_phasespace)


ifig=0;




ifig=ifig+1;
figure(ifig)
hist(beam_phasespace(1,:),20)
xlabel('X [m]','FontSize',16)
filename = ['fig' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);


ifig=ifig+1;
figure(ifig)
hist(beam_phasespace(2,:),20)
xlabel('X prime','FontSize',16)
filename = ['fig' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);

ifig=ifig+1;
figure(ifig)
hist(beam_phasespace(3,:),20)
xlabel('Y [m]','FontSize',16)
filename = ['fig' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);


ifig=ifig+1;
figure(ifig)
hist(beam_phasespace(4,:),20)
xlabel('Y prime','FontSize',16)
filename = ['fig' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);

ifig=ifig+1;
figure(ifig)
hist(beam_phasespace(5,:),20)
xlabel('T [s]','FontSize',16)
filename = ['fig' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);


ifig=ifig+1;
figure(ifig)
hist(beam_phasespace(6,:),20)
xlabel('delta','FontSize',16)
filename = ['fig' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);