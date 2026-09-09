clc; clear all;

if ~exist('Stab1D','dir'); mkdir('Stab1D'); end

Nx = 128;
Lx = 0; Rx = 1; h = (Rx-Lx)/Nx;
x = linspace(Lx+0.5*h,Rx-0.5*h,Nx)';

m = 12; eps = m*h/(2*sqrt(2)*atanh(0.9));

dt = 0.5*h^2*eps^2/(4*h^2+2*eps^2);
Nt = 100;
T  = Nt*dt;

phi_ave = 0;
alpha0  = 0.001;

Kstar = 1/(pi*eps);
Klist = 1:ceil(1.6*Kstar);

sigma_ana = 1/eps^2 - (Klist*pi).^2;
sigma_num = zeros(size(Klist));

for k = 1:numel(Klist)
    K = Klist(k);
    phi0 = phi_ave + alpha0*cos(K*pi*x);
    phi = [phi0(1); phi0; phi0(end)];

    for it = 1:Nt
        phi(1) = phi(2); phi(Nx+2) = phi(Nx+1);
        lam = sum(phi(2:Nx+1)-phi(2:Nx+1).^3,'all')/sum(1-phi(2:Nx+1).^2,'all');
        phi(2:Nx+1) = phi(2:Nx+1)+dt*((phi(2:Nx+1)-phi(2:Nx+1).^3)/eps^2 ...
            +(phi(1:Nx)+phi(3:Nx+2)-2*phi(2:Nx+1))/h^2-lam*(1-phi(2:Nx+1).^2)/eps^2);
    end

    amp = max(abs(phi(2:Nx+1)-phi_ave),[],'all');
    sigma_num(k) = log(amp/alpha0)/T;
end

figure(1); clf; hold on; box on; grid on;
set(gca,'fontsize',21);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'PaperPositionMode','auto','Position',[500 500 800 400])
plot(Klist,sigma_ana,'k-','linewidth',2)
plot(Klist,sigma_num,'ko','markersize',10,'LineWidth',1.5)
leg = legend('Analytic $\sigma(K)$','Numerical $\tilde\sigma(K)$');
set(leg,'Interpreter','latex','Location','northeast')
text('Interpreter','latex','string','$K$','FontSize',23,'Position',[21.8409645027594 -4394.43078690269])
text('Interpreter','latex','string','$\sigma$','FontSize',27,'Position',[-0.474975850872115 1115.12117156829])
axis([1 23 -4000 2200])
drawnow
print('Stab1D/Fig2.eps','-depsc','-loose')