clc; clear all; close all;

if ~exist('Convergence','dir'); mkdir('Convergence'); end

mlist = 1:6;
T_fac = 2;
ratio = 1;
m_ref = 12;

%% 1D
Nx = 2048;
Lx = 0; Rx = 1; h = (Rx-Lx)/Nx;
x  = linspace(Lx+0.5*h,Rx-0.5*h,Nx)';

m   = 12;
eps = m*h/(2*sqrt(2)*atanh(0.9));

dt_max = h^2*eps^2/(4*h^2+2*eps^2);
T = T_fac*dt_max;

dtList = dt_max./2.^mlist';
NtList = round(T./dtList);

Nx_ref = ratio*Nx; h_ref = (Rx-Lx)/Nx_ref;
x_ref  = linspace(Lx+0.5*h_ref,Rx-0.5*h_ref,Nx_ref)';
dt_ref = dt_max/2^m_ref;
Nt_ref = round(T/dt_ref);

phi0_ref = 0.1*cos(pi*x_ref);
phiRef = solve_CAC_1D(phi0_ref,Nx_ref,h_ref,eps,dt_ref,Nt_ref);

phi0 = 0.1*cos(pi*x);
Error = zeros(numel(mlist),1);
for l = 1:numel(mlist)
    phi = solve_CAC_1D(phi0,Nx,h,eps,dtList(l),NtList(l));
    Error(l) = sqrt(sum((phiRef-phi).^2,'all')/Nx);
end
Order = log2(Error(1:end-1)./Error(2:end));

fprintf('\n1D convergence test, Nx = %d, dt_ref = %.6e, T = %.6e\n',Nx,dt_ref,T);
fprintf('%-16s %-16s %-10s\n','dt','Error','Order');
for l = 1:numel(mlist)
    if l == 1
        fprintf('%-16.6e %-16.4e %-10s\n',dtList(l),Error(l),'-');
    else
        fprintf('%-16.6e %-16.4e %-10.2f\n',dtList(l),Error(l),Order(l-1));
    end
end

dtList1D = dtList;
Error1D = Error;

figure(1); clf; hold on; box on; grid on;
set(gca,'fontsize',21);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'PaperPositionMode','auto','Position',[500 500 550 400])
loglog(dtList1D,Error1D,'ko-','markersize',10,'linewidth',1.5,'markerfacecolor','k')
loglog([dtList1D(4) dtList1D(3) dtList1D(3) dtList1D(4)],[Error1D(4)*0.65 2*(Error1D(4)*0.65) Error1D(4)*0.65 Error1D(4)*0.65] ,'k-','markersize',10,'linewidth',1.5,'markerfacecolor','k')
set(gca,'XScale','log','YScale','log')
text('Interpreter','latex','String','$\Vert {\bf e}^{\Delta t} \Vert_2$','Position',[5.05581402385649e-10 3.5707401120964e-05],'FontSize',18)
text('Interpreter','latex','String','$\Delta t$','Position',[4.42111334172424e-08 1.35106223314442e-06],'FontSize',18)
text('Interpreter','latex','String','$1$','Position',[1.27587642296974e-08 9.91758380326531e-06],'FontSize',18)
text('Interpreter','latex','String','$1$','Position',[8.10837014915596e-09 5.66936615801347e-06],'FontSize',18)
axis([0.65*dtList1D(end) 1.4*dtList1D(1) 0.65*Error1D(end) 1.4*Error1D(1)])
xticks([1.0e-9 1.0e-8 1.0e-7 1.0e-6])
drawnow
print('Convergence/convergence_1D.eps','-depsc','-loose')

%% 2D
Nx = 512; Ny = 512;
Lx = 0; Rx = 1; h = (Rx-Lx)/Nx;
Ly = 0; Ry = 1;
x  = linspace(Lx+0.5*h,Rx-0.5*h,Nx)';
y  = linspace(Ly+0.5*h,Ry-0.5*h,Ny)';

m   = 12;
eps = m*h/(2*sqrt(2)*atanh(0.9));

dt_max = h^2*eps^2/(4*h^2+4*eps^2);
T = T_fac*dt_max;

dtList = dt_max./2.^mlist';
NtList = round(T./dtList);

Nx_ref = ratio*Nx; Ny_ref = ratio*Ny; h_ref = (Rx-Lx)/Nx_ref;
x_ref  = linspace(Lx+0.5*h_ref,Rx-0.5*h_ref,Nx_ref)';
y_ref  = linspace(Ly+0.5*h_ref,Ry-0.5*h_ref,Ny_ref)';
dt_ref = dt_max/2^m_ref;
Nt_ref = round(T/dt_ref);

phi0_ref = 0.1*(cos(pi*x_ref)*cos(pi*y_ref)');
phiRef = solve_CAC_2D(phi0_ref,Nx_ref,Ny_ref,h_ref,eps,dt_ref,Nt_ref);

phi0 = 0.1*(cos(pi*x)*cos(pi*y)');
Error = zeros(numel(mlist),1);
for l = 1:numel(mlist)
    phi = solve_CAC_2D(phi0,Nx,Ny,h,eps,dtList(l),NtList(l));
    Error(l) = sqrt(sum((phiRef-phi).^2,'all')/(Nx*Ny));
end
Order = log2(Error(1:end-1)./Error(2:end));

fprintf('\n2D convergence test, Nx = Ny = %d, dt_ref = %.6e, T = %.6e\n',Nx,dt_ref,T);
fprintf('%-16s %-16s %-10s\n','dt','Error','Order');
for l = 1:numel(mlist)
    if l == 1
        fprintf('%-16.6e %-16.4e %-10s\n',dtList(l),Error(l),'-');
    else
        fprintf('%-16.6e %-16.4e %-10.2f\n',dtList(l),Error(l),Order(l-1));
    end
end

dtList2D = dtList;
Error2D = Error;

figure(2); clf; hold on; box on; grid on;
set(gca,'fontsize',21);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'PaperPositionMode','auto','Position',[500 500 550 400])
loglog(dtList2D,Error2D,'ko-','markersize',10,'linewidth',1.5,'markerfacecolor','k')
E0 = Error2D(4)*0.65;
loglog([dtList2D(4) dtList2D(3) dtList2D(3) dtList2D(4)],[E0 2*E0 E0 E0],'k-','markersize',10,'linewidth',1.5,'markerfacecolor','k')
set(gca,'XScale','log','YScale','log')
text('Interpreter','latex','String','$\Vert {\bf e}^{\Delta t} \Vert_2$','Position',[4.48698892492482e-09 2.02325276491976e-05],'FontSize',18)
text('Interpreter','latex','String','$\Delta t$','Position',[4.22027340858541e-07 2.72669934596015e-07],'FontSize',18)
text('Interpreter','latex','String','$1$','Position',[7.1631258401484e-08 1.19639702220196e-06],'FontSize',18)
text('Interpreter','latex','String','$1$','Position',[1.13082073192321e-07 2.05992570322413e-06],'FontSize',18)
axis([0.65*dtList2D(end) 1.4*dtList2D(1) 0.65*Error2D(end) 1.4*Error2D(1)])
xticks([1.0e-9 1.0e-8 1.0e-7 1.0e-6])
drawnow
print('Convergence/convergence_2D.eps','-depsc','-loose')

%% 3D
Nx = 128; Ny = 128; Nz = 128;
Lx = 0; Rx = 1; h = (Rx-Lx)/Nx;
Ly = 0; Ry = 1;
Lz = 0; Rz = 1;
x  = linspace(Lx+0.5*h,Rx-0.5*h,Nx)';
y  = linspace(Ly+0.5*h,Ry-0.5*h,Ny)';
z  = linspace(Lz+0.5*h,Rz-0.5*h,Nz)';

m   = 12;
eps = m*h/(2*sqrt(2)*atanh(0.9));

dt_max = h^2*eps^2/(4*h^2+6*eps^2);
T = T_fac*dt_max;

dtList = dt_max./2.^mlist';
NtList = round(T./dtList);

Nx_ref = ratio*Nx; Ny_ref = ratio*Ny; Nz_ref = ratio*Nz; h_ref = (Rx-Lx)/Nx_ref;
x_ref  = linspace(Lx+0.5*h_ref,Rx-0.5*h_ref,Nx_ref)';
y_ref  = linspace(Ly+0.5*h_ref,Ry-0.5*h_ref,Ny_ref)';
z_ref  = linspace(Lz+0.5*h_ref,Rz-0.5*h_ref,Nz_ref)';
dt_ref = dt_max/2^m_ref;
Nt_ref = round(T/dt_ref);

phi0_ref = 0.1*(cos(pi*x_ref) .* reshape(cos(pi*y_ref),1,Ny_ref,1) .* reshape(cos(pi*z_ref),1,1,Nz_ref));
phiRef = solve_CAC_3D(phi0_ref,Nx_ref,Ny_ref,Nz_ref,h_ref,eps,dt_ref,Nt_ref);

phi0 = 0.1*(cos(pi*x) .* reshape(cos(pi*y),1,Ny,1) .* reshape(cos(pi*z),1,1,Nz));
Error = zeros(numel(mlist),1);
for l = 1:numel(mlist)
    phi = solve_CAC_3D(phi0,Nx,Ny,Nz,h,eps,dtList(l),NtList(l));
    Error(l) = sqrt(sum((phiRef-phi).^2,'all')/(Nx*Ny*Nz));
end
Order = log2(Error(1:end-1)./Error(2:end));

fprintf('\n3D convergence test, Nx = Ny = Nz = %d, dt_ref = %.6e, T = %.6e\n',Nx,dt_ref,T);
fprintf('%-16s %-16s %-10s\n','dt','Error','Order');
for l = 1:numel(mlist)
    if l == 1
        fprintf('%-16.6e %-16.4e %-10s\n',dtList(l),Error(l),'-');
    else
        fprintf('%-16.6e %-16.4e %-10.2f\n',dtList(l),Error(l),Order(l-1));
    end
end

dtList3D = dtList;
Error3D = Error;

figure(3); clf; hold on; box on; grid on;
set(gca,'fontsize',21);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'PaperPositionMode','auto','Position',[500 500 550 400])
loglog(dtList3D,Error3D,'ko-','markersize',10,'linewidth',1.5,'markerfacecolor','k')
E0 = 0.65*Error3D(4);
loglog([dtList3D(4) dtList3D(3) dtList3D(3) dtList3D(4)],[E0 2*E0 E0 E0],'k-','markersize',10,'linewidth',1.5,'markerfacecolor','k')
set(gca,'XScale','log','YScale','log')
text('Interpreter','latex','String','$\Vert {\bf e}^{\Delta t} \Vert_2$','Position',[4.98097738076717e-08 6.51999761684217e-06],'FontSize',18)
text('Interpreter','latex','String','$\Delta t$','Position',[4.24226120063779e-06 9.20932819346808e-08],'FontSize',18)
text('Interpreter','latex','String','$1$','Position',[7.99935004980046e-07 3.9049278123423e-07],'FontSize',18)
text('Interpreter','latex','String','$1$','Position',[1.25093033964779e-06 6.7181796775103e-07],'FontSize',18)
axis([0.65*dtList3D(end) 1.4*dtList3D(1) 0.65*Error3D(end) 1.4*Error3D(1)])
xticks([1.0e-9 1.0e-8 1.0e-7 1.0e-6 1.0e-5])
drawnow
print('Convergence/convergence_3D.eps','-depsc','-loose')

function phi_T = solve_CAC_1D(phi0,Nx,h,eps,dt,Nt)

phi = [phi0(1); phi0; phi0(end)];

for it = 1:Nt
phi(1) = phi(2); phi(Nx+2) = phi(Nx+1);

lam = sum(phi(2:Nx+1)-phi(2:Nx+1).^3,'all')/sum(1-phi(2:Nx+1).^2,'all');

phi(2:Nx+1) = phi(2:Nx+1)+dt*((phi(2:Nx+1)-phi(2:Nx+1).^3)/eps^2 ...
    +(phi(1:Nx)+phi(3:Nx+2)-2*phi(2:Nx+1))/h^2 ...
    -lam*(1-phi(2:Nx+1).^2)/eps^2);
end
phi_T = phi(2:Nx+1);
end

function phi_T = solve_CAC_2D(phi0,Nx,Ny,h,eps,dt,Nt)
phi = zeros(Nx+2,Ny+2);
phi(2:Nx+1,2:Ny+1) = phi0;

for it = 1:Nt
phi(1,:) = phi(2,:); phi(Nx+2,:) = phi(Nx+1,:);
phi(:,1) = phi(:,2); phi(:,Ny+2) = phi(:,Ny+1);

lam = sum(phi(2:Nx+1,2:Ny+1)-phi(2:Nx+1,2:Ny+1).^3,'all')/sum(1-phi(2:Nx+1,2:Ny+1).^2,'all');
phi(2:Nx+1,2:Ny+1) = phi(2:Nx+1,2:Ny+1)+dt*((phi(2:Nx+1,2:Ny+1)-phi(2:Nx+1,2:Ny+1).^3)/eps^2 ...
    +(phi(1:Nx,2:Ny+1)+phi(3:Nx+2,2:Ny+1)+phi(2:Nx+1,1:Ny)+phi(2:Nx+1,3:Ny+2) ...
    -4*phi(2:Nx+1,2:Ny+1))/h^2 ...
    -lam*(1-phi(2:Nx+1,2:Ny+1).^2)/eps^2);
end
phi_T = phi(2:Nx+1,2:Ny+1);
end

function phi_T = solve_CAC_3D(phi0,Nx,Ny,Nz,h,eps,dt,Nt)
phi = zeros(Nx+2,Ny+2,Nz+2);
phi(2:Nx+1,2:Ny+1,2:Nz+1) = phi0;

for it = 1:Nt
phi(1,:,:) = phi(2,:,:); phi(Nx+2,:,:) = phi(Nx+1,:,:);
phi(:,1,:) = phi(:,2,:); phi(:,Ny+2,:) = phi(:,Ny+1,:);
phi(:,:,1) = phi(:,:,2); phi(:,:,Nz+2) = phi(:,:,Nz+1);

lam = sum(phi(2:Nx+1,2:Ny+1,2:Nz+1)-phi(2:Nx+1,2:Ny+1,2:Nz+1).^3,'all')/sum(1-phi(2:Nx+1,2:Ny+1,2:Nz+1).^2,'all');
phi(2:Nx+1,2:Ny+1,2:Nz+1) = phi(2:Nx+1,2:Ny+1,2:Nz+1)+dt*((phi(2:Nx+1,2:Ny+1,2:Nz+1)-phi(2:Nx+1,2:Ny+1,2:Nz+1).^3)/eps^2 ...
    +(phi(1:Nx,2:Ny+1,2:Nz+1)+phi(3:Nx+2,2:Ny+1,2:Nz+1) ...
    +phi(2:Nx+1,1:Ny,2:Nz+1)+phi(2:Nx+1,3:Ny+2,2:Nz+1) ...
    +phi(2:Nx+1,2:Ny+1,1:Nz)+phi(2:Nx+1,2:Ny+1,3:Nz+2) ...
    -6*phi(2:Nx+1,2:Ny+1,2:Nz+1))/h^2 ...
    -lam*(1-phi(2:Nx+1,2:Ny+1,2:Nz+1).^2)/eps^2);
end
phi_T = phi(2:Nx+1,2:Ny+1,2:Nz+1);
end
