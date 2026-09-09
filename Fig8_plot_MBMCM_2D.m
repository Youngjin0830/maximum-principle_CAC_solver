clc; clear; close all;

tex2 = sprintf('MBMCM2/MBMCM2_%dit.mat',0);
load(tex2);

figure(1); clf; hold on; box on;
set(gca,'fontsize',21);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto')
contour(x,y,phi',[0 0],'k','LineWidth',1)
axis image
axis([Lx Rx Ly Ry])
drawnow

figure(2); clf; hold on; box on;
set(gca,'fontsize',21);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto')
surf(x(2:Nx+1),y(2:Ny+1),phi(2:Nx+1,2:Ny+1)','EdgeColor','none')
clim([-1 1]);
view(2)
colormap jet
shading interp
axis image
axis([Lx Rx Ly Ry])
text('Interpreter','latex','String','$x$','FontSize',24,'Position',[0.751337726149568 -0.0608218069243655])
text('Interpreter','latex','String','$y$','FontSize',24,'Position',[-0.0755105228776694 0.917777414865518])
drawnow
tex = sprintf('MBMCM2/MBMCM2_%dit.eps',0);
print(tex,'-depsc')

lin = {'--','-.',':'};
count = 0;
for it = [600 1400 5000]
    count = count+1;
    tex2 = sprintf('MBMCM2/MBMCM2_%dit.mat',it);
    load(tex2,'phi');

    figure(it); clf; hold on; box on;
    set(gca,'fontsize',21);
    set(gca, 'TickLabelInterpreter', 'latex');
    set(gcf,'PaperPositionMode','auto')
    surf(x(2:Nx+1),y(2:Ny+1),phi(2:Nx+1,2:Ny+1)','EdgeColor','none')
    clim([-1 1]);
    view(2)
    colormap jet
    shading interp
    axis image
    axis([Lx Rx Ly Ry])
    text('Interpreter','latex','String','$x$','FontSize',24,'Position',[0.751337726149568 -0.0608218069243655])
    text('Interpreter','latex','String','$y$','FontSize',24,'Position',[-0.0755105228776694 0.917777414865518])
    drawnow
    tex = sprintf('MBMCM2/MBMCM2_%dit.eps',it);
    print(tex,'-depsc')

    figure(1); 
    set(gca,'fontsize',21);
    set(gca, 'TickLabelInterpreter', 'latex');
    set(gcf,'PaperPositionMode','auto')
    contour(x,y,phi',[0 0],'k','LineStyle',lin{count},'LineWidth',1)
    axis image
    axis([Lx Rx Ly Ry])
    drawnow
end
text('Interpreter','latex','String','$x$','FontSize',24,'Position',[0.751337726149568 -0.0608218069243655])
text('Interpreter','latex','String','$y$','FontSize',24,'Position',[-0.0755105228776694 0.917777414865518])

leg = legend('$t=0$','$t=600\Delta t$','$t=1400\Delta t$','$t=5000\Delta t$');
set(leg,'interpreter','latex')

annotation('arrow',[0.362204724409449 0.428696412948381],...
    [0.168253968253968 0.315873015873016])
annotation('arrow',[0.483814523184602 0.388451443569554],...
    [0.544444444444445 0.634920634920635])

tex = sprintf('MBMCM2/MBMCM2_cont.eps');
print(tex,'-depsc')
