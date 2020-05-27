function [Am,T]= p_matrix3D(grid,par,mu)
%
% Computes flow system matrix.
%
  Nx = grid.Nx;
  Ny = grid.Ny;
  Nz = grid.Nz;
  dx = grid.dx;
  dy = grid.dy;
  dz = grid.dz;
  dx2 = dx*dx;
  dy2 = dy*dy;
  dz2 = dz*dz;
%  
  Np = Nx*Ny*Nz; 
%
% Transmisibility matrices.
% Tx = 2(dy*dz)/(1/K1 +1/K2)
% so that
% q = -Tx*(p1-p0)/dx
%
  T.Tx = zeros(Nx+1,Ny,Nz);
  T.Tx(2:Nx,:,:) = (2*dy*dz)./(mu(1:Nx-1,:,:) + mu(2:Nx,:,:));
%
  T.Ty = zeros(Nx,Ny+1,Nz);

  T.Ty(:,2:Ny,:) = (2*dx*dz)./(mu(:,1:Ny-1,:) + mu(:,2:Ny,:));
%
  T.Tz = zeros(Nx,Ny,Nz+1);
  T.Tz(:,:,2:Nz) = (2*dx*dy)./(mu(:,:,1:Nz-1) + mu(:,:,2:Nz));

%
  T.Tx1 = reshape(T.Tx(1:Nx,:,:),Np,1);
  T.Tx2 = reshape(T.Tx(2:Nx+1,:,:),Np,1);
%
  T.Ty1 = reshape(T.Ty(:,1:Ny,:),Np,1);
  T.Ty2 = reshape(T.Ty(:,2:Ny+1,:),Np,1);
%
  T.Tz1 = reshape(T.Tz(:,:,1:Nz),Np,1);
  T.Tz2 = reshape(T.Tz(:,:,2:Nz+1),Np,1);
%
%
% Assemble system of equations
%
  Dp = zeros(Np,1);
  Dp = T.Tx1/dx2 + T.Ty1/dy2 + T.Tz1/dz2 + T.Tx2/dx2 + T.Ty2/dy2 + T.Tz2/dz2;

%
%  if par.LeftIsDirichletFlow
 
      T.TzDirichL = reshape((2*dx*dy).*(1./mu(:,:,1)),Nx*Ny,1);

      Dp(1:Nx*Ny,1) = ...
                   T.Tx1(1:Nx*Ny,1)/dx2 + ...                 
                   T.Ty1(1:Nx*Ny,1)/dy2 + ...
                   T.Tx2(1:Nx*Ny,1)/dx2 + ...
                   T.Ty2(1:Nx*Ny,1)/dy2 + ...
                   T.Tz2(1:Nx*Ny,1)/dz2 + ...
                   T.TzDirichL/dz2;
%
%  if par.RightIsDirichletFlow
 
      T.TzDirichR = reshape((2*dx*dy).*(1./mu(:,:,Nz)),Nx*Ny,1);
      
      Dp(Np-(Nx*Ny)+1:Np,1) = ... 
                         T.Tx1(Np-(Nx*Ny)+1:Np,1)/dx2 + ...
                         T.Ty1(Np-(Nx*Ny)+1:Np,1)/dy2 + ...
                         T.Tx2(Np-(Nx*Ny)+1:Np,1)/dx2 + ...
                         T.Ty2(Np-(Nx*Ny)+1:Np,1)/dy2 + ...
                         T.Tz1(Np-(Nx*Ny)+1:Np,1)/dz2 + ...
                         T.TzDirichR/dz2;

%

    Am = spdiags(...
     [-T.Tz2/dz2, ...
      -T.Ty2/dy2, ...
      -T.Tx2/dx2,  Dp, -T.Tx1/dx2,  ...
      -T.Ty1/dy2, ...
      -T.Tz1/dz2], ...
   [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny], Np, Np);
% 
end
