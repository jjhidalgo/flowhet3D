close all;
clear all;
% heterogeneity
var_lnk= 1;%5.9;
corr_lenx= 1./5;
corr_leny= 1./5;
corr_lenz= 1./5;
Nz = 32;%2^9;
%
% Set to true to plot logK, pressure, etc. If false only
% histograms are plot.
plotit = true; 
%
% Geometry (Lz=1 always).
%
grid.A = 1;     %A = Lx/Lz
grid.W = 1;     %W = Ly/Lz
%
% Discretization
%
grid.Nz = Nz;
grid.Nx = round(grid.A*grid.Nz);
grid.Ny = round(grid.W*grid.Nz);

grid.Lx = grid.A;
grid.Ly = grid.W;
grid.Lz = 1;
grid.dx = grid.Lx/grid.Nx;
grid.dy = grid.Ly/grid.Ny;
grid.dz = grid.Lz/grid.Nz;
%
dx = grid.dx;
dy = grid.dy;
dz = grid.dz;
Lx = grid.Lx;
Ly = grid.Ly;
Lz = grid.Lz;
Nx = grid.Nx;
Ny = grid.Ny;
Nz = grid.Nz;

%
% Boundary condtions
%
% z=0 & z=Lz are Dirichlet boundaries.
par.pL = 1.;
par.pR = 0.;


% modes
kx = (2*pi/Lx)*[0:(Nx/2-1) (-Nx/2):(-1)]';
ky = (2*pi/Ly)*[0:(Ny/2-1) (-Ny/2):(-1)]';
kz = (2*pi/Lz)*[0:(Nz/2-1) (-Nz/2):(-1)]';
[kx,ky,kz]= meshgrid(kx,ky,kz);

[K.kperm,var_lnk_actual, mean_lnk_actual]= gen_randperm3D(var_lnk, ...
                                   corr_lenx,corr_leny,corr_lenz,kx,ky,kz);
                               
disp(strcat(['Actual mean-logK = ' num2str(mean_lnk_actual)]));
disp(strcat(['Actual var-logK = ' num2str(var_lnk_actual)]));


if plotit
  figure(1)
  h=slice(log(K.kperm),1:Nx,1:Ny,1:Nz);
  set(h, 'EdgeColor','none', 'FaceColor','interp');
  axis equal
  title('log-K')
end


% Solves flow
[Am,Trans] = p_matrix3D(grid,par,1./K.kperm);
[ux,uz,uy,p] = p_rhs3D(grid,par,Am,Trans);

if plotit
  figure(2)
  h=slice(p,1:Nx,1:Ny,1:Nz);
  set(h, 'EdgeColor','none', 'FaceColor','interp');
  colorbar()
  title('p')
end

vel = sqrt(ux.*ux + uy.*uy + uz.*uz);


figure(3)
hv = histogram(log(vel(:)), 'Normalization', 'pdf');
hold on
hk = histogram(log(K.kperm(:)), 'Normalization', 'pdf');
title('histograms log-v; log-K')

figure(5)
edges = hv.BinEdges;
edgesv = (edges(1:end-1)+edges(2:end))/2.;
plot(edgesv, hv.Values, 'b-', 'DisplayName', 'vel')
hold on
edges = hk.BinEdges;
edgesk = (edges(1:end-1)+edges(2:end))/2.;
plot(edgesk, hk.Values, 'r-', 'DisplayName', 'k')
%set(gca,'xscale', 'log')
%set(gca,'yscale', 'log')
legend('location', 'best')
title('histograms log-v; log-K')
%saves data
dlmwrite('histV.dat', [edgesv' hv.Values'],'delimiter', ' ');
dlmwrite('histK.dat', [edgesk' hk.Values'],'delimiter', ' ');

return
