clc; clear all; close all;

if ~exist('EnergyDiss1D','dir'); mkdir('EnergyDiss1D'); end

rng(7)

Nx = 256;
Lx = 0; Rx = 1; h = (Rx-Lx)/Nx;
x = linspace(Lx-0.5*h,Rx+0.5*h,Nx+2);

phi0 = 0.1*(2*rand(Nx+2,1)-1);
phi = phi0;

m = 8; eps = m*h/(2*sqrt(2)*atanh(0.9));
dt = 1*h^2*eps^2/(4*h^2+2*eps^2);

T = 0.04;
Nt = ceil(T/dt);

energy = zeros(Nt+1,1);
E1D = @(phi_sol,Nx,h,eps) ...
    (1/(2*h))*sum((phi_sol(2:Nx+2)-phi_sol(1:Nx+1)).^2,'all') ...
    +(h/(4*eps^2))*sum((1-phi_sol(2:Nx+1).^2).^2,'all');
energy(1) = E1D(phi,Nx,h,eps);

for it = 1:Nt
    phi(1) = phi(2); phi(Nx+2) = phi(Nx+1);

    f = phi(2:Nx+1)-phi(2:Nx+1).^3;
    C = sum(f,'all')/sum((1-phi(2:Nx+1).^2),'all');
    phi(2:Nx+1) = phi(2:Nx+1)+dt*(f/eps^2 ...
        +(phi(1:Nx)+phi(3:Nx+2)-2*phi(2:Nx+1))/h^2-C*(1-phi(2:Nx+1).^2)/eps^2);

    energy(it+1) = E1D(phi,Nx,h,eps);
end
t = [0:Nt]*dt;
energy_normalized_log = energy/energy(2);

figure(2); clf;
semilogx(t,energy_normalized_log,'b-','linewidth',1.5); hold on; box on; grid on;
set(gca,'fontsize',21);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto','Position',[500 500 500 420])
text('Interpreter','latex','String','$t$','FontSize',24,'Units','normalized','Position',[0.961671145778262 -0.074591439688716])
text('Interpreter','latex','String','$\frac{\mathcal{E}_h(\phi^n)}{\mathcal{E}_h(\phi^1)}$','FontSize',17.5,'Units','normalized','Position',[-0.181122809808032 0.777548638132296])
xticks([1.0e-5 1.0e-4 1.0e-3 1.0e-2 1.0e-1])
axis([t(2) 0.04 0 1])
yticks([0 0.5 1])
drawnow
print('EnergyDiss1D/normalized_energy_1D_log.eps','-depsc')
