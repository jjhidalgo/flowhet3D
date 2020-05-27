function [perm,var_lnk_actual, mean_lnk_actual]=gen_randperm3D(var_lnk,lx,ly,lz,kx,ky,kz)
i= sqrt(-1);
[Nx,Ny,Nz] = size(kx);
dkx = kx(1,2,1);
dky = ky(2,1,1);
dkz = kz(1,1,2);
%
%uncomment next 1 line to obtain unique random numbers everytime matlab is started
%RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

%%%%Spectral density functions%%%%
%Whittle-A spectrum, isotropic medium
%a=pi/(4*sqrt(corr_lenx*corr_leny));
%S=(2*var_lnk*a*a/pi)*( (KX2+KY2) ./ ((a*a+KX2+KY2).^3) );  

%Gaussian spectrum, Random Field Generator, Ruan and McLaughlin, Adv Water Res. 1998,
%S=(1/2/pi)*var_lnk*corr_lenx*corr_leny*exp(-0.5*(corr_lenx^2)*(KX2) - 0.5*(corr_leny^2)*(KY2));

%Spectrum of modified exp autocovariance, Gelhar & Axness, WRR 1983
S = (var_lnk*lx*ly*lz)./(pi*pi*(1 + (lx^2).*kx.^2 + (ly^2).*ky.^2 + (lz^2).*kz.^2).^2); 

H=S.^(0.5); 
theta = 2*pi*rand(Nx,Ny,Nz); %random phase angle

dZ = H.*exp(i*theta).*sqrt(dkx*dky*dkz)*Nx*Ny*Nz;
lnK=ifftn(dZ); %log(hyd. cond.)
K1=sqrt(2)*real(lnK); K2=sqrt(2)*imag(lnK); 
var_lnk_actual=var(reshape(K1,Nx*Ny*Nz,1));
mean_lnk_actual= mean(reshape(K1,Nx*Ny*Nz,1));
K1=exp(K1); K2=exp(K2);

%choose either K1 or K2
perm = K1;

end
