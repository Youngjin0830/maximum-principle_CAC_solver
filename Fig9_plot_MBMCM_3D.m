clc; clear; close all;

fs = 26;

tex2 = sprintf('MBMCM3/MBMCM3_%dit.mat',0);
load(tex2);
figure(1); clf; hold on; box on; grid on;
set(gca,'fontsize',fs);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto')
plotphi = permute(phi,[2 1 3]);
s1 = isosurface(xx,yy,zz,plotphi,0);
patch(s1,'facecolor',		"#FF0000",'facealpha',1,'edgecolor','none');
axis image
view([-10 20]);
camlight('right'); daspect([1 1 1]); lighting phong;
yticks([0 1])
axis([Lx Rx Ly Ry Lz Rz])
text('Interpreter','latex','String','$x$','Position',[1.17310183064304 -0.306643996176831 -0.0171531714917474],'FontSize',fs+2)
text('Interpreter','latex','String','$y$','Position',[-0.0882982711809737 0.747433456247258 -0.124639614324941],'FontSize',fs)
text('Interpreter','latex','String','$z$','Position',[-0.166264055177213 0.787216242393029 0.908201030882068],'FontSize',fs+2)
drawnow
tex = sprintf('MBMCM3/MBMCM3_%dit.eps',0);
print(tex,'-depsc')

for it = [300 500 1200]
    tex2 = sprintf('MBMCM3/MBMCM3_%dit.mat',it);
    load(tex2,'phi');

    figure(it); clf; hold on; box on; grid on;
    set(gca,'fontsize',fs);
    set(gca, 'TickLabelInterpreter', 'latex');
    set(gcf,'PaperPositionMode','auto')
    plotphi = permute(phi,[2 1 3]);
    s1 = isosurface(xx,yy,zz,plotphi,0);
    patch(s1,'facecolor',		"#FF0000",'facealpha',1,'edgecolor','none');
    axis image
    view([-10 20]);
    camlight('right'); daspect([1 1 1]); lighting phong;
    yticks([0 1])
    axis([Lx Rx Ly Ry Lz Rz])
    text('Interpreter','latex','String','$x$','Position',[1.17310183064304 -0.306643996176831 -0.0171531714917474],'FontSize',fs+2)
    text('Interpreter','latex','String','$y$','Position',[-0.0882982711809737 0.747433456247258 -0.124639614324941],'FontSize',fs)
    text('Interpreter','latex','String','$z$','Position',[-0.166264055177213 0.787216242393029 0.908201030882068],'FontSize',fs+2)
    drawnow
    tex = sprintf('MBMCM3/MBMCM3_%dit.eps',it);
    print(tex,'-depsc')
end

it = 1200;
tex2 = sprintf('MBMCM3/MBMCM3_max_min_mass.mat_%dit.mat',it);
load(tex2);
t = [0:it]*dt;
ms = 9;
figure(999); clf; hold on; box on; grid on;
set(gca,'fontsize',fs);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto')
plot(t,maxi,'k-^','markersize',ms,'markerfacecolor','k','linewidth',1,'MarkerIndices',[1:100:it+1])
plot(t,mini,'k-v','markersize',ms,'markerfacecolor','k','linewidth',1,'MarkerIndices',[1:100:it+1])
plot(t,mass,'k-o','markersize',ms,'markerfacecolor','k','linewidth',1,'MarkerIndices',[1:100:it+1])
axis([0 t(1201) -1.2 1.2])
xticks([0 0.002:0.002:0.6])
xticklabels({'$0$','$0.002$','$0.004$','$0.006$'})
text('Interpreter','latex','String','$t$','FontSize',28,'Position',[0.003 -1.37647058823529])
leg = legend('$Max(\phi^n)$','$Min(\phi^n)$','$\mathcal{M}(\phi^n)$');
set(leg,'Interpreter','latex','Location','east','Position',[0.632080647422181 0.23931954920977 0.180316888396157 0.269396551724138])
print('MBMCM3/maximum_test_3D.eps','-depsc')
