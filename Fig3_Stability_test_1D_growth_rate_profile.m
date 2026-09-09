clc; clear all;

if ~exist('Stab1D','dir'); mkdir('Stab1D'); end

Nx = 128;
Lx = 0; Rx = 1; h = (Rx-Lx)/Nx;
x = linspace(Lx+0.5*h,Rx-0.5*h,Nx)';

m = 12; eps = m*h/(2*sqrt(2)*atanh(0.9));
dt = 0.5*h^2*eps^2/(4*h^2+2*eps^2);

alpha0 = 0.001;
K = 4;
sigma_K = 1/eps^2-(K*pi)^2;

Nt = 100;

phi0 = alpha0*cos(K*pi*x);
phi = [phi0(1); phi0; phi0(end)];

u_linear = alpha0*cos(K*pi*x);

figure(1); clf; hold on; box on; grid on;
set(gca,'fontsize',21,'TickLabelInterpreter','latex');
set(gcf,'PaperPositionMode','auto','Position',[500 500 550 500])
plot(x,u_linear,'k-','linewidth',2);
plot(x,phi(2:Nx+1),'k*','linewidth',1.5,'markersize',15,'MarkerIndices',1:6:Nx);
axis([0 1 -1e-2 1e-2])
text('Interpreter','latex','string','$x$','FontSize',25,'Position',[0.739059967585089 -0.0109983633387889])
text('Interpreter','latex','string','$\phi$','FontSize',25,'Position',[-0.109400324149109 0.00765957446808508])
leg = legend('Linear $u(x,t)$', 'Numerical $\phi(x,t)$');
set(leg, 'Interpreter', 'latex', 'Location', 'south', 'FontSize', 20,'Position',[0.440133522727273 0.146 0.443712121212121 0.137]);
tex = sprintf('Stab1D/Fig3_%dit.eps',0);
print(tex,'-depsc')

for it = 1:Nt
phi(1) = phi(2); phi(Nx+2) = phi(Nx+1);

lam = sum(phi(2:Nx+1)-phi(2:Nx+1).^3,'all')/sum(1-phi(2:Nx+1).^2,'all');

phi(2:Nx+1) = phi(2:Nx+1)+dt*((phi(2:Nx+1)-phi(2:Nx+1).^3)/eps^2 ...
    +(phi(1:Nx)+phi(3:Nx+2)-2*phi(2:Nx+1))/h^2-lam*(1-phi(2:Nx+1).^2)/eps^2);

if it == 60 || it == 100
    t = it*dt;
    u_linear = alpha0*exp(sigma_K*t)*cos(K*pi*x);

    figure; clf; hold on; box on; grid on;
    set(gca,'fontsize',21,'TickLabelInterpreter','latex');
    set(gcf,'PaperPositionMode','auto','Position',[500 500 550 500])
    plot(x,u_linear,'k-','linewidth',2);
    plot(x,phi(2:Nx+1),'k*','linewidth',1.5,'markersize',15,'MarkerIndices',1:6:Nx);
    axis([0 1 -1e-2 1e-2])
    text('Interpreter','latex','string','$x$','FontSize',25,'Position',[0.739059967585089 -0.0109983633387889])
    text('Interpreter','latex','string','$\phi$','FontSize',25,'Position',[-0.109400324149109 0.00765957446808508])
    leg = legend('Linear $u(x,t)$', 'Numerical $\phi(x,t)$');
    set(leg, 'Interpreter', 'latex', 'Location', 'south', 'FontSize', 20,'Position',[0.440133522727273 0.146 0.443712121212121 0.137]);

    tex = sprintf('Stab1D/Fig3_%dit.eps',it);
    print(tex,'-depsc')
end
end