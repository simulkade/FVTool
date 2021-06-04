function faceGrad = gradientTermCylindrical3D(phi)
% this function calculates the gradient of a variable in x and y direction
% it checks for the availability of the ghost variables and use them, otherwise
% estimate them, assuming a zero gradient on the boundaries
%
% SYNOPSIS:
%   faceGrad = gradientTermCylindrical3D(phi)
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

% check the size of the variable and the mesh dimension
Nr = phi.domain.dims(1);
Ntheta = phi.domain.dims(2);
Nz = phi.domain.dims(3);
DR = repmat(phi.domain.cellsize.x, 1, Ntheta, Nz);
DTHETA = repmat(phi.domain.cellsize.y', Nr, 1, Nz);
DZ = zeros(1,1,Nz+2);
DZ(1,1,:) = phi.domain.cellsize.z;
DZ = repmat(DZ, Nr, Ntheta, 1);
dr = 0.5*(DR(1:end-1,:,:)+DR(2:end,:,:));
dtetta = 0.5*(DTHETA(:,1:end-1,:)+DTHETA(:,2:end,:));
dz = 0.5*(DZ(:,:,1:end-1)+DZ(:,:,2:end));
rp = repmat(phi.domain.cellcenters.x, 1, Ntheta+1, Nz);

xvalue = (phi.value(2:Nr+2,2:Ntheta+1,2:Nz+1)-phi.value(1:Nr+1,2:Ntheta+1,2:Nz+1))./dr;
yvalue = (phi.value(2:Nr+1,2:Ntheta+2,2:Nz+1)-phi.value(2:Nr+1,1:Ntheta+1,2:Nz+1))./(dtetta.*rp);
zvalue = (phi.value(2:Nr+1,2:Ntheta+1,2:Nz+2)-phi.value(2:Nr+1,2:Ntheta+1,1:Nz+1))./dz;

faceGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
