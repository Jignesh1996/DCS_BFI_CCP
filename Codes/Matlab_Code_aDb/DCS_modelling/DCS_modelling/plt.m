function plt(x,y,fname,xlbl,ylbl,lgnd)

hfig = figure()
plot(x,y);
fname = fname;
legend(lgnd);

% legend("Normalized 1cm","Norm 2.5cm","Norm 2.5cm TR","DCS 1cm","DCS 2.5cm","TR sys DCS 2.5cm ")
title(fname)
xlabel(xlbl);
ylabel(ylbl);
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-dpng','-painters')




end