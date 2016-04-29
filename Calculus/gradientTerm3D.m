function faceGrad = gradientTerm3D(phi)
% this function calculates the gradient of a variable in x and y direction
% it checks for the availability of the ghost variables and use them, otherwise
% estimate them, assuming a zero gradient on the boundaries
%
% SYNOPSIS:
%
%
% PARAMETERS:
%
%
% RETURNS:
%
%
% EXAMPLE:
%
% SEE ALSO:
%

% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file

% check the size of the variable and the mesh dimension
Nx = phi.domain.dims(1);
Ny = phi.domain.dims(2);
Nz = phi.domain.dims(3);
DX = repmat(phi.domain.cellsize.x, 1, Ny, Nz);
DY = repmat(phi.domain.cellsize.y', Nx, 1, Nz);
DZ = zeros(1,1,Nz+2);
DZ(1,1,:) = phi.domain.cellsize.z;
DZ=repmat(DZ, Nx, Ny, 1);
dx = 0.5*(DX(1:end-1,:,:)+DX(2:end,:,:));
dy = 0.5*(DY(:,1:end-1,:)+DY(:,2:end,:));
dz = 0.5*(DZ(:,:,1:end-1)+DZ(:,:,2:end));


% in this case, ghost cells have values
xvalue = (phi.value(2:Nx+2,2:Ny+1,2:Nz+1)-phi.value(1:Nx+1,2:Ny+1,2:Nz+1))./dx;
yvalue = (phi.value(2:Nx+1,2:Ny+2,2:Nz+1)-phi.value(2:Nx+1,1:Ny+1,2:Nz+1))./dy;
zvalue = (phi.value(2:Nx+1,2:Ny+1,2:Nz+2)-phi.value(2:Nx+1,2:Ny+1,1:Nz+1))./dz;


faceGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
