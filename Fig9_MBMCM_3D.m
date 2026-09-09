clc; clear; close all;

if ~exist('MBMCM3','dir'); mkdir('MBMCM3'); end

Nx = 200*1.5;
Lx = 0; Rx = 1.5; h = (Rx-Lx)/Nx; 
Ly = 0; Ry = 1; Ny = (Ry-Ly)/h;
Lz = 0; Rz = 1; Nz = (Rz-Lz)/h;

x = linspace(Lx-0.5*h,Rx+0.5*h,Nx+2);
y = linspace(Ly-0.5*h,Ry+0.5*h,Ny+2);
z = linspace(Lz-0.5*h,Rz+0.5*h,Nz+2);
[xx,yy,zz]=meshgrid(x,y,z);

m = 8;
eps = m*h/(2*sqrt(2)*atanh(0.9));

X_shifted = 1*(xx-0.75); Y_shifted = 1*(yy-0.5); Z_shifted = 1*(zz-0.5);
R = sqrt(X_shifted.^2+Y_shifted.^2+Z_shifted.^2);
Theta = acos(X_shifted./(R+1.0e-19));
Phi = zeros(size(xx));

idx1 = (Y_shifted > 0); idx2 = ~idx1;

Phi(idx1) = atan(Z_shifted(idx1)./Y_shifted(idx1));
Phi(idx2) = pi+atan(Z_shifted(idx2)./Y_shifted(idx2));
Phi(isnan(Phi)) = 0;

l = 16; m = 4;
P_lm = zeros(size(Theta));

for i = 1:numel(Theta)
    P_all = legendre(l,cos(Theta(i)),'norm');
    P_lm(i) = P_all(m+1);
end

Y_16_4 = P_lm.*cos(m * Phi);

phi = tanh((0.3+0.2*Y_16_4-R)/(sqrt(2)*eps));
phi = permute(phi,[2 1 3]);

axisfs = 19; fs=25; lw=1;
figure(1); clf; hold on; box on;
set(gca,'fontsize',21);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto')
plotphi = permute(phi,[2 1 3]);
s1 = isosurface(xx,yy,zz,plotphi,0);
patch(s1,'facecolor',	"#0072BD",'facealpha',1,'edgecolor','none');
axis image
view([-30 15]);
camlight('right'); daspect([1 1 1]); lighting phong;
axis([Lx Rx Ly Ry Lz Rz])
text('Interpreter','latex','String','$x$','Position',[1.14279740663677 -0.296518174715338 -0.00421347836988417],'FontSize',fs)
text('Interpreter','latex','String','$y$','Position',[-0.148122476038889 0.796117106303228 -0.0214561373730593],'FontSize',fs-2)
text('Interpreter','latex','String','$z$','Position',[-0.274032821775638 0.791310573359389 0.867863173010198],'FontSize',fs)
drawnow
tex = sprintf('MBMCM3/MBMCM3_%dit.eps',0);
print(tex,'-depsc')

dt = 1*h^2*eps^2/(4*h^2+6*eps^2);
tex2 = sprintf('MBMCM3/MBMCM3_%dit.mat',0);
save(tex2);
Nt = 1200;

mass = zeros(Nt,1);
maxi = mass;
mini = mass;
mass(1) = sum((phi(2:Nx+1,2:Ny+1,2:Nz+1)+1)/2,"all")*h^3;
maxi(1) = max(phi,[],'all');
mini(1) = min(phi,[],'all');

for it = 1:Nt

phi(1,:,:) = phi(2,:,:); phi(Nx+2,:,:) = phi(Nx+1,:,:);
phi(:,1,:) = phi(:,2,:); phi(:,Ny+2,:) = phi(:,Ny+1,:);
phi(:,:,1) = phi(:,:,2); phi(:,:,Nz+2) = phi(:,:,Nz+1);

lam = sum((phi(2:Nx+1,2:Ny+1,2:Nz+1)-phi(2:Nx+1,2:Ny+1,2:Nz+1).^3),'all') ...
    /sum((1-phi(2:Nx+1,2:Ny+1,2:Nz+1).^2),'all');

phi(2:Nx+1,2:Ny+1,2:Nz+1) = phi(2:Nx+1,2:Ny+1,2:Nz+1)+dt*( ...
    (phi(2:Nx+1,2:Ny+1,2:Nz+1)-phi(2:Nx+1,2:Ny+1,2:Nz+1).^3)/eps^2 ...
    +(phi(1:Nx,2:Ny+1,2:Nz+1)+phi(3:Nx+2,2:Ny+1,2:Nz+1) ...
    +phi(2:Nx+1,1:Ny,2:Nz+1)+phi(2:Nx+1,3:Ny+2,2:Nz+1) ...
    +phi(2:Nx+1,2:Ny+1,1:Nz)+phi(2:Nx+1,2:Ny+1,3:Nz+2) ...
    -6*phi(2:Nx+1,2:Ny+1,2:Nz+1))/h^2 ...
    -lam.*(1-phi(2:Nx+1,2:Ny+1,2:Nz+1).^2)/eps^2);

if it == 300 || it == 500 || it == 1200
    figure(1); clf; hold on; box on;
    set(gca,'fontsize',21);
    set(gca, 'TickLabelInterpreter', 'latex');
    set(gcf,'PaperPositionMode','auto')
    plotphi = permute(phi,[2 1 3]);
    s1 = isosurface(xx,yy,zz,plotphi,0);
    patch(s1,'facecolor',	"#0072BD",'facealpha',1,'edgecolor','none');
    axis image
    view([-30 15]);
    camlight('right'); daspect([1 1 1]); lighting phong;
    axis([Lx Rx Ly Ry Lz Rz])
    text('Interpreter','latex','String','$x$','Position',[1.14279740663677 -0.296518174715338 -0.00421347836988417],'FontSize',fs)
    text('Interpreter','latex','String','$y$','Position',[-0.148122476038889 0.796117106303228 -0.0214561373730593],'FontSize',fs-2)
    text('Interpreter','latex','String','$z$','Position',[-0.274032821775638 0.791310573359389 0.867863173010198],'FontSize',fs)
    drawnow
    tex2 = sprintf('MBMCM3/MBMCM3_%dit.mat',it);
    save(tex2,'phi');
end
maxi(it+1) = max(phi,[],'all');
mini(it+1) = min(phi,[],'all');
mass(it+1) = sum((phi(2:Nx+1,2:Ny+1,2:Nz+1)+1)/2,"all")*h^3;
end
tex2 = sprintf('MBMCM3/MBMCM3_max_min_mass.mat_%dit.mat',it);
save(tex2,'maxi','mini','mass');
