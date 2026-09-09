clc; clear all; close all;

if ~exist('EnergyDiss3D','dir'); mkdir('EnergyDiss3D'); end

rng(7)

Nx = 256;
Lx = 0; Rx = 1; h = (Rx-Lx)/Nx;
Ly = 0; Ry = 1; Ny = (Ry-Ly)/h;
Lz = 0; Rz = 1; Nz = (Rz-Lz)/h;

x = linspace(Lx-0.5*h,Rx+0.5*h,Nx+2);
y = linspace(Ly-0.5*h,Ry+0.5*h,Ny+2);
z = linspace(Lz-0.5*h,Rz+0.5*h,Nz+2);

phi0 = 0.1*(2*rand(Nx+2,Ny+2,Nz+2)-1);
phi = phi0;

m = 8; eps = m*h/(2*sqrt(2)*atanh(0.9));
dt = 1*h^2*eps^2/(4*h^2+6*eps^2);

T = 0.04;
Nt = ceil(T/dt);

energy = zeros(Nt+1,1);
E3D = @(phi_sol,Nx,Ny,Nz,h,eps) ...
    (h/2)*(sum((phi_sol(2:Nx+2,2:Ny+1,2:Nz+1)-phi_sol(1:Nx+1,2:Ny+1,2:Nz+1)).^2,'all') ...
    +sum((phi_sol(2:Nx+1,2:Ny+2,2:Nz+1)-phi_sol(2:Nx+1,1:Ny+1,2:Nz+1)).^2,'all') ...
    +sum((phi_sol(2:Nx+1,2:Ny+1,2:Nz+2)-phi_sol(2:Nx+1,2:Ny+1,1:Nz+1)).^2,'all')) ...
    +(h^3/(4*eps^2))*sum((1-phi_sol(2:Nx+1,2:Ny+1,2:Nz+1).^2).^2,'all');
energy(1) = E3D(phi,Nx,Ny,Nz,h,eps);

for it = 1:Nt
    phi(1,:,:) = phi(2,:,:); phi(Nx+2,:,:) = phi(Nx+1,:,:);
    phi(:,1,:) = phi(:,2,:); phi(:,Ny+2,:) = phi(:,Ny+1,:);
    phi(:,:,1) = phi(:,:,2); phi(:,:,Nz+2) = phi(:,:,Nz+1);

    f = phi(2:Nx+1,2:Ny+1,2:Nz+1)-phi(2:Nx+1,2:Ny+1,2:Nz+1).^3;
    C = sum(f,'all')/sum((1-phi(2:Nx+1,2:Ny+1,2:Nz+1).^2),'all');
    lap = (phi(1:Nx,2:Ny+1,2:Nz+1)+phi(3:Nx+2,2:Ny+1,2:Nz+1) ...
        +phi(2:Nx+1,1:Ny,2:Nz+1)+phi(2:Nx+1,3:Ny+2,2:Nz+1) ...
        +phi(2:Nx+1,2:Ny+1,1:Nz)+phi(2:Nx+1,2:Ny+1,3:Nz+2) ...
        -6*phi(2:Nx+1,2:Ny+1,2:Nz+1))/h^2;
    phi(2:Nx+1,2:Ny+1,2:Nz+1) = phi(2:Nx+1,2:Ny+1,2:Nz+1)+dt*(f/eps^2+lap ...
        -C*(1-phi(2:Nx+1,2:Ny+1,2:Nz+1).^2)/eps^2);

    energy(it+1) = E3D(phi,Nx,Ny,Nz,h,eps);
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
print('EnergyDiss3D/normalized_energy_3D_log.eps','-depsc')
