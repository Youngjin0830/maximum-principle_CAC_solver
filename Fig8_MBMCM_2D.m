clc; clear; close all;

if ~exist('MBMCM2','dir'); mkdir('MBMCM2'); end

Nx = 256;
Lx = 0; Rx = 1; h = (Rx-Lx)/Nx; 
Ly = 0; Ry = 1; Ny = (Ry-Ly)/h;

x = linspace(Lx-0.5*h,Rx+0.5*h,Nx+2);
y = linspace(Ly-0.5*h,Ry+0.5*h,Ny+2);
[X, Y] = meshgrid(x, y);
m = 11;
eps = m*h/(2*sqrt(2)*atanh(0.9));

R = 0.3;
per = 0.15;
n_theta = 7;
xc = 0.5; yc = 0.5;
theta = atan2(Y-yc,X-xc);
r = sqrt((X-xc).^2+(Y-yc).^2);
R_perturbed = R+per*sin(n_theta*theta);
phi = -tanh((r-R_perturbed)/(sqrt(2)*eps));

figure(1); clf; hold on; box on;
set(gca,'fontsize',21);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto')
contourf(x,y,phi,[0 0])
colormap jet
shading interp
axis image
axis([Lx Rx Ly Ry])
text('Interpreter','latex','String','$x$','FontSize',24,'Position',[0.742582862336339 -0.084168110426311])
text('Interpreter','latex','String','$y$','FontSize',24,'Position',[-0.103599221789883 1.28550583657588])
drawnow

dt = 1*h^2*eps^2/(4*h^2+4*eps^2);
tex2 = sprintf('MBMCM2/MBMCM2_%dit.mat',0);
save(tex2);
Nt = 5000;

for it = 1:Nt
phi(1,:) = phi(2,:); phi(Nx+2,:) = phi(Nx+1,:);
phi(:,1) = phi(:,2); phi(:,Ny+2) = phi(:,Ny+1);

lam = sum((phi(2:Nx+1,2:Ny+1)-phi(2:Nx+1,2:Ny+1).^3),'all')/sum((1-phi(2:Nx+1,2:Ny+1).^2),'all');
phi(2:Nx+1,2:Ny+1) = phi(2:Nx+1,2:Ny+1)+dt*((phi(2:Nx+1,2:Ny+1)-phi(2:Nx+1,2:Ny+1).^3)/eps^2 ...
    +(phi(1:Nx,2:Ny+1)+phi(3:Nx+2,2:Ny+1)+phi(2:Nx+1,1:Ny)+phi(2:Nx+1,3:Ny+2)-4*phi(2:Nx+1,2:Ny+1))/h^2 ...
    -lam.*(1-phi(2:Nx+1,2:Ny+1).^2)/eps^2);

if it == 600 || it == 1400 || it == 5000
    figure(1); clf; hold on; box on;
    set(gca,'fontsize',21);
    set(gca, 'TickLabelInterpreter', 'latex');
    set(gcf,'PaperPositionMode','auto')
    contourf(x,y,phi,[0 0])
    colormap jet
    shading interp
    axis image
    axis([Lx Rx Ly Ry])
    text('Interpreter','latex','String','$x$','FontSize',24,'Position',[0.742582862336339 -0.084168110426311])
    text('Interpreter','latex','String','$y$','FontSize',24,'Position',[-0.103599221789883 1.28550583657588])
    drawnow
    tex2 = sprintf('MBMCM2/MBMCM2_%dit.mat',it);
    save(tex2,'phi');
end
end
