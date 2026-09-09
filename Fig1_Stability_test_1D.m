clc; clear;

if ~exist('Stab1D','dir'); mkdir('Stab1D'); end

rng(7)

Nx = 128; Lx = 0; Rx = 1; h = (Rx-Lx)/Nx; 
x = linspace(Lx-0.5*h,Rx+0.5*h,Nx+2);
phi0 = 0.1*(2*rand(Nx+2,1)-1);
m = 12; eps = m*h/(2*sqrt(2)*atanh(0.9));

%% 1.2\Delta t_{max}
dt = 1.2*h^2*eps^2/(4*h^2+2*eps^2);
phi = phi0;
Nt = 1500;

mass = zeros(Nt,1); maxi = mass; mini = mass;
mass(1) = sum((phi(2:Nx+1)+1)/2,"all")*h;
maxi(1) = max(phi,[],'all');
mini(1) = min(phi,[],'all');

for it = 1:Nt

phi(1) = phi(2); phi(Nx+2) = phi(Nx+1);
phi(2:Nx+1) = phi(2:Nx+1)+dt*((phi(2:Nx+1)-phi(2:Nx+1).^3)/eps^2 ...
    +(phi(1:Nx)+phi(3:Nx+2)-2*phi(2:Nx+1))/h^2 ...
    -sum((phi(2:Nx+1)-phi(2:Nx+1).^3),'all')/sum((1-phi(2:Nx+1).^2),'all').*(1-phi(2:Nx+1).^2)/eps^2);

if it == 100 || it == 500 || it == 1500
    figure(1); clf; hold on; box on; grid on;
    set(gca,'fontsize',21);
    set(gca, 'TickLabelInterpreter', 'latex');
    set(gcf,'PaperPositionMode','auto','Position',[500 500 450 400])
    plot(x,phi,'linewidth',1)
    axis([Lx Rx -1.5 1.5])
    text('Interpreter','latex','String','$x$','FontSize',24,'Position',[0.735269128368905 -1.69358303335679])
    text('Interpreter','latex','String','$\phi$','FontSize',24,'Position',[-0.0900492085472616 1.27587809290842])
    drawnow
    tex = sprintf('Stab1D/Fig1_b_%dit.eps',it);
    print(tex,'-depsc')
end
maxi(it+1) = max(phi,[],'all');
mini(it+1) = min(phi,[],'all');
mass(it+1) = sum((phi(2:Nx+1)+1)/2,"all")*h;
end
maxi2 = maxi; mini2 = mini; mass2 = mass;
t2 = [0:it]*dt; Nt2 = Nt;

%% \Delta t_{max}
dt = 1*h^2*eps^2/(4*h^2+2*eps^2);
phi = phi0;
Nt = Nt*1.2;

mass = zeros(Nt,1); maxi = mass; mini = mass;
mass(1) = sum((phi(2:Nx+1)+1)/2,"all")*h;
maxi(1) = max(phi,[],'all');
mini(1) = min(phi,[],'all');

for it = 1:Nt

phi(1) = phi(2); phi(Nx+2) = phi(Nx+1);
phi(2:Nx+1) = phi(2:Nx+1)+dt*((phi(2:Nx+1)-phi(2:Nx+1).^3)/eps^2 ...
    +(phi(1:Nx)+phi(3:Nx+2)-2*phi(2:Nx+1))/h^2 ...
    -sum((phi(2:Nx+1)-phi(2:Nx+1).^3),'all')/sum((1-phi(2:Nx+1).^2),'all').*(1-phi(2:Nx+1).^2)/eps^2);

if it == 120 || it == 600 || it == 1800
    figure(2); clf; hold on; box on; grid on;
    set(gca,'fontsize',21);
    set(gca, 'TickLabelInterpreter', 'latex');
    set(gcf,'PaperPositionMode','auto','Position',[500 500 450 400])
    plot(x,phi,'linewidth',1)
    axis([Lx Rx -1.5 1.5])
    text('Interpreter','latex','String','$x$','FontSize',24,'Position',[0.735269128368905 -1.69358303335679])
    text('Interpreter','latex','String','$\phi$','FontSize',24,'Position',[-0.0900492085472616 1.27587809290842])
    drawnow
    tex = sprintf('Stab1D/Fig1_a_%dit.eps',it);
    print(tex,'-depsc')
end
maxi(it+1) = max(phi,[],'all');
mini(it+1) = min(phi,[],'all');
mass(it+1) = sum((phi(2:Nx+1)+1)/2,"all")*h;
end
t = [0:it]*dt;

ms = 9;
figure(3); clf; hold on; box on; grid on;
set(gca,'fontsize',21);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto')
plot(t,maxi,'b-^','markersize',ms,'markerfacecolor','b','linewidth',1,'MarkerIndices',[1:100:Nt])
plot(t,mini,'b-v','markersize',ms,'markerfacecolor','b','linewidth',1,'MarkerIndices',[1:100:Nt])
plot(t,mass,'b-o','markersize',ms,'markerfacecolor','b','linewidth',1,'MarkerIndices',[1:100:Nt])
plot(t2,maxi2,'r--^','markersize',ms,'markerfacecolor','r','linewidth',1,'MarkerIndices',[50:100:Nt2])
plot(t2,mini2,'r--v','markersize',ms,'markerfacecolor','r','linewidth',1,'MarkerIndices',[50:100:Nt2])
plot(t2,mass2,'r--o','markersize',ms,'markerfacecolor','r','linewidth',1,'MarkerIndices',[50:100:Nt2])
axis([0 t(end) -1.4 1.4])
text('Interpreter','latex','String','$t$','FontSize',24,'Position',[0.0428089532613149 -1.5692269629563])
leg = legend('$Max(\phi^n),~\Delta t = \Delta t_{max}$','$Min(\phi^n),~\Delta t = \Delta t_{max}$','$\mathcal{M}(\phi^n),~\Delta t = \Delta t_{max}$','$Max(\phi^n),~\Delta t = 1.2\Delta t_{max}$','$Min(\phi^n),~\Delta t = 1.2\Delta t_{max}$','$\mathcal{M}(\phi^n),~\Delta t = 1.2\Delta t_{max}$');
set(leg,'Interpreter','latex','Location','east')
print('Stab1D/Fig1_c.eps','-depsc')
