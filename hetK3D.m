close all;
clear all;
% heterogeneity
var_lnk = 1;%5.9;
lx = 1./5;
ly = 1./5;
lz = 1./5;
N = 2^8;
L = 1.;
dx = L/N;

%Check discretization
Nmin = 5*(1/lx); %5 cells per correlation length.
straux = strcat(['N = ' num2str(N) '. Minimum Nx = ' num2str(Nmin)]);
if Nmin>N
     straux = strcat([straux ' (not ok)']);
else
    straux = strcat([straux ' (ok)']);
end
disp(straux)
%


% modes
k = (2*pi/L)*[0:(N/2-1) (-N/2):(-1)]';
[kx,ky,kz] = meshgrid(k,k,k);

[K.kperm,var_lnk_actual]= gen_randperm3D(var_lnk,lx,ly,lz,kx,ky,kz);
disp(strcat(['Actual var-log = ' num2str(var_lnk_actual)]));

figure(1)
h=slice(log(K.kperm),1:N,1:N,1:N);
set(h, 'EdgeColor','none', 'FaceColor','interp');
title('log-K')
return
