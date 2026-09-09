clc; clear all; close all;

if ~exist('Stab2D','dir'); mkdir('Stab2D'); end

rng(123)

Nx = 256; Lx = 0; Rx = 1; h = (Rx-Lx)/Nx; 
Ly = 0; Ry = 1; Ny = (Ry-Ly)/h;
x = linspace(Lx-0.5*h,Rx+0.5*h,Nx+2);
y = linspace(Ly-0.5*h,Ry+0.5*h,Ny+2);
phi = 0.5*(2*rand(Nx+2,Ny+2)-1);
m = 8; eps = m*h/(2*sqrt(2)*atanh(0.9));
dt = 1*h^2*eps^2/(4*h^2+4*eps^2);
Nt = 1000;

mass = zeros(Nt,1); maxi = mass; mini = mass;
mass(1) = sum((phi(2:Nx+1,2:Ny+1)+1)/2,"all")*h^2;
maxi(1) = max(phi,[],'all');
mini(1) = min(phi,[],'all');

for it = 1:Nt

phi(1,:) = phi(2,:); phi(Nx+2,:) = phi(Nx+1,:);
phi(:,1) = phi(:,2); phi(:,Ny+2) = phi(:,Ny+1);
phi(2:Nx+1,2:Ny+1) = phi(2:Nx+1,2:Ny+1)+dt*((phi(2:Nx+1,2:Ny+1)-phi(2:Nx+1,2:Ny+1).^3)/eps^2+(phi(1:Nx,2:Ny+1)+phi(3:Nx+2,2:Ny+1)+phi(2:Nx+1,1:Ny)+phi(2:Nx+1,3:Ny+2)-4*phi(2:Nx+1,2:Ny+1))/h^2 ...
    -sum((phi(2:Nx+1,2:Ny+1)-phi(2:Nx+1,2:Ny+1).^3),'all')/sum((1-phi(2:Nx+1,2:Ny+1).^2),'all').*(1-phi(2:Nx+1,2:Ny+1).^2)/eps^2);

if it == 50 || it == 200 || it == 1000
    figure(1); clf; hold on; box on;
    set(gca,'fontsize',21);
    set(gca, 'TickLabelInterpreter', 'latex');
    set(gcf,'PaperPositionMode','auto')
    surf(x(2:Nx+1),y(2:Ny+1),phi(2:Nx+1,2:Ny+1),'EdgeColor','none')
    clim([-1 1]); colormap jet; shading interp; axis image;
    axis([Lx Rx Ly Ry])
    text('Interpreter','latex','String','$x$','FontSize',24,'Position',[0.744528387628164 -0.05838990030958])
    text('Interpreter','latex','String','$y$','FontSize',24,'Position',[-0.0729912806384733 0.919699508648357])
    drawnow
    tex = sprintf('Stab2D/Fig4a_%dit.eps',it);
    print(tex,'-depsc')
end
maxi(it+1) = max(phi,[],'all');
mini(it+1) = min(phi,[],'all');
mass(it+1) = sum((phi(2:Nx+1,2:Ny+1)+1)/2,"all")*h^2;
end
t = [0:it]*dt;

%%
ms = 9;
figure(3); clf; hold on; box on; grid on;
set(gca,'fontsize',21);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto')
plot(t,maxi,'b-^','markersize',ms,'markerfacecolor','b','linewidth',1,'MarkerIndices',[1:100:Nt+1])
plot(t,mini,'b-v','markersize',ms,'markerfacecolor','b','linewidth',1,'MarkerIndices',[1:100:Nt+1])
plot(t,mass,'b-o','markersize',ms,'markerfacecolor','b','linewidth',1,'MarkerIndices',[1:100:Nt+1])
axis([0 t(end) -1.01 1.01])
text('Interpreter','latex','String','$t$','FontSize',24,'Position',[0.003 -0])
xticks([0 0.001 0.002 0.003])
xticklabels({'$0$','$0.001$','$0.002$','$0.003$'})
leg = legend('$Max(\phi^n)$','$Min(\phi^n)$','$\mathcal{M}(\phi^n)$');
set(leg,'Interpreter','latex','Location','east')
print('Stab2D/Fig4b.eps','-depsc')
