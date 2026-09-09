clc; clear all; close all;

if ~exist('MBM1','dir'); mkdir('MBM1'); end

Nx = 128;
Lx = 0; Rx = 1; h = (Rx-Lx)/Nx; 
Ly = 0; Ry = 1.5; Ny = (Ry-Ly)/h;

x = linspace(Lx-0.5*h,Rx+0.5*h,Nx+2);
y = linspace(Ly-0.5*h,Ry+0.5*h,Ny+2);
m = 8; eps = m*h/(2*sqrt(2)*atanh(0.9));

x1 = 0.5; y1 = 1.15;
x2 = 0.5; y2 = 0.5;
phi = 1+tanh((0.2-sqrt((x'-x1).^2+(y-y1).^2))/(sqrt(2)*eps))+tanh((0.3-sqrt((x'-x2).^2+(y-y2).^2))/(sqrt(2)*eps));

xc = [x2; x1]; yc = [y2; y1];
it = 0;
c = contour(x,y,phi',[0 0],'EdgeColor','none');
k = 1; count = 1;
while k-1 ~= size(c,2)
    k2 = k+c(2,k);
    radius(it+1,count) = mean(sqrt(sum((c(:,k+1:k2)-[xc(count); yc(count)]).^2,1)));
    count = count+1;
    k = k2+1;
end

figure(1); clf; hold on; box on;
set(gca,'fontsize',21);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto')
surf(x(2:Nx+1),y(2:Ny+1),phi(2:Nx+1,2:Ny+1)','EdgeColor','none')
clim([-1 1]);
colormap jet
shading interp
axis image
axis([Lx Rx Ly Ry])
text('Interpreter','latex','String','$x$','FontSize',24,'Position',[0.742582862336339 -0.084168110426311])
text('Interpreter','latex','String','$y$','FontSize',24,'Position',[-0.103599221789883 1.28550583657588])
drawnow
tex = sprintf('MBM1/MBMM1_%dit.eps',it);
print(tex,'-depsc')

dt = 1*h^2*eps^2/(4*h^2+4*eps^2);
Nt = 5500;

for it = 1:Nt
phi(1,:) = phi(2,:); phi(Nx+2,:) = phi(Nx+1,:);
phi(:,1) = phi(:,2); phi(:,Ny+2) = phi(:,Ny+1);

lam = sum((phi(2:Nx+1,2:Ny+1)-phi(2:Nx+1,2:Ny+1).^3),'all')/sum((1-phi(2:Nx+1,2:Ny+1).^2),'all');
phi(2:Nx+1,2:Ny+1) = phi(2:Nx+1,2:Ny+1)+dt*((phi(2:Nx+1,2:Ny+1)-phi(2:Nx+1,2:Ny+1).^3)/eps^2 ...
    +(phi(1:Nx,2:Ny+1)+phi(3:Nx+2,2:Ny+1)+phi(2:Nx+1,1:Ny)+phi(2:Nx+1,3:Ny+2)-4*phi(2:Nx+1,2:Ny+1))/h^2 ...
    -lam.*(1-phi(2:Nx+1,2:Ny+1).^2)/eps^2);

c = contour(x,y,phi',[0 0],'EdgeColor','none');
k = 1; count = 1;
while k-1 ~= size(c,2)
    k2 = k+c(2,k);
    radius(it+1,count) = mean(sqrt(sum((c(:,k+1:k2)-[xc(count); yc(count)]).^2,1)));
    count = count+1;
    k = k2+1;
end

if it == 3000 || it == 4000 || it == 4500
    figure(1); clf; hold on; box on;
    set(gca,'fontsize',21);
    set(gca, 'TickLabelInterpreter', 'latex');
    set(gcf,'PaperPositionMode','auto')
    surf(x(2:Nx+1),y(2:Ny+1),phi(2:Nx+1,2:Ny+1)','EdgeColor','none')
    clim([-1 1]);
    colormap jet
    shading interp
    axis image
    axis([Lx Rx Ly Ry])
    text('Interpreter','latex','String','$x$','FontSize',24,'Position',[0.742582862336339 -0.084168110426311])
    text('Interpreter','latex','String','$y$','FontSize',24,'Position',[-0.103599221789883 1.28550583657588])
    drawnow
    tex = sprintf('MBM1/MBMM1_%dit.eps',it);
    print(tex,'-depsc')
end
end

r0 = 0.2; R0 = 0.3;
tf = 0.055; N = 5500;
r(1) = r0; R(1) = R0;

for n = 1:N
    rn = r(n);
    Rn = R(n);

    k1r = 2 / (rn + Rn) - 1 / rn;
    k1R = 2 / (rn + Rn) - 1 / Rn;

    r_half = rn + 0.5 * dt * k1r;
    R_half = Rn + 0.5 * dt * k1R;
    k2r = 2 / (r_half + R_half) - 1 / r_half;
    k2R = 2 / (r_half + R_half) - 1 / R_half;

    r_half = rn + 0.5 * dt * k2r;
    R_half = Rn + 0.5 * dt * k2R;
    k3r = 2 / (r_half + R_half) - 1 / r_half;
    k3R = 2 / (r_half + R_half) - 1 / R_half;

    r_full = rn + dt * k3r;
    R_full = Rn + dt * k3R;
    k4r = 2 / (r_full + R_full) - 1 / r_full;
    k4R = 2 / (r_full + R_full) - 1 / R_full;

    r(n+1) = rn + (dt / 6) * (k1r + 2*k2r + 2*k3r + k4r);
    R(n+1) = Rn + (dt / 6) * (k1R + 2*k2R + 2*k3R + k4R);
    if r(n+1) <0 || R(n+1) <0
        break;
    end
end
r(n+2:N+1) = r(n+1);
R(n+2:N+1) = R(n+1);

lw = 1; ms = 8;
t = [0:it]*dt; t2 = [0:N]*dt;
figure(3); clf; hold on; box on; grid on;
set(gca,'fontsize',17);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto')
plot(t2,R,'k--','linewidth',lw)
plot(t2,r,'k-','linewidth',lw)
plot(t,radius(:,1),'ro','linewidth',lw,'markersize',ms,'MarkerIndices',[1:300:4401 4450])
plot(t,radius(:,2),'b*','linewidth',lw,'markersize',ms,'MarkerIndices',[1:300:4401 4450])
leg = legend('$R_{ref}$','$r_{ref}$','$R^n$','$r^n$');
set(leg,'interpreter','latex')
text('Interpreter','latex','String','$t$','FontSize',19,'Position',[0.0449490373725934 -0.0321513002364066])
axis([0 0.0535 0 0.4])
print('radius.eps','-depsc')
