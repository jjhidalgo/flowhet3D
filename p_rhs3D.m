function [ux_avg,uy_avg,uz_avg,p]= p_rhs3D(grid,par,Am,T)

% Computes right hand side of flow system and solves.
%
%

  Nx = grid.Nx;
  Ny = grid.Ny;
  Nz = grid.Nz;
  Np = Nx*Ny*Nz;
  dx = grid.dx;
  dy = grid.dy;
  dz = grid.dz;
  dx2 = dx*dx;
  dy2 = dy*dy;
  dz2 = dz*dz;
  
 
%
% Source term
%
%
  S = zeros(Np,1);

%
% Boundary conditions
%
  S(1:Nx*Ny,1) = S(1:Nx*Ny,1) + T.TzDirichL.*(par.pL/dz2);
  
  S(Np-(Nx*Ny)+1:Np,1) = S(Np-(Nx*Ny)+1:Np,1) + T.TzDirichR.*(par.pR/dz2);

%
% Solve system of equations
%
  p = Am\S;
  p = reshape(p,Nx,Ny,Nz);
  p = full(p);
  
  %Am = full(Am);f
  %keyboard()
% Compute velocities
%
  ux = zeros(Nx+1,Ny,Nz);
  uy = zeros(Nx,Ny+1,Nz);
  uz = zeros(Nx,Ny,Nz+1);

%
%

  ux(2:Nx,:,:) =  -T.Tx(2:Nx,:,:).*(p(2:Nx,:,:)-p(1:Nx-1,:,:))/dx;
  uy(:,2:Ny,:) =  -T.Ty(:,2:Ny,:).*(p(:,2:Ny,:)-p(:,1:Ny-1,:))/dy;
  
  uz(:,:,2:Nz) = -T.Tz(:,:,2:Nz).*(p(:,:,2:Nz) - p(:,:,1:Nz-1))/dz;
  

  %if par.LeftIsDirichletFlow
  uz(:,:,1) = -reshape(T.TzDirichL,Nx,Ny).*(p(:,:,1)-par.pL)/(dx);
    
  
  %if par.RightIsDirichletFlow
  uz(:,:,Nz+1) = -reshape(T.TzDirichR,Nx,Ny).*(par.pR-p(:,:,Nz))/(dx);

%
% Average velocity in the cell.
%
  
  ux_avg = 0.5*( ux(1:Nx,:,:) + ux(2:Nx+1,:,:) )/(dy*dz);
  uy_avg = 0.5*( uy(:,1:Ny,:) + uy(:,2:Ny+1,:) )/(dx*dz);
  uz_avg = 0.5*( uz(:,:,1:Nz) + uz(:,:,2:Nz+1) )/(dx*dy);
%
end
