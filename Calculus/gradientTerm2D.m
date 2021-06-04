function faceGrad = gradientTerm2D(phi)
% this function calculates the gradient of a variable in x and y direction
% it checks for the availability of the ghost variables and use them, otherwise
% estimate them, assuming a zero gradient on the boundaries
%
% SYNOPSIS:
%   faceGrad = gradientTerm2D(phi)
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
Nx = phi.domain.dims(1);
Ny = phi.domain.dims(2);
DX = repmat(phi.domain.cellsize.x, 1, Ny);
DY = repmat(phi.domain.cellsize.y', Nx, 1);
dx = 0.5*(DX(1:end-1,:)+DX(2:end,:));
dy = 0.5*(DY(:,1:end-1)+DY(:,2:end));
% in this case, ghost cells have values
xvalue = (phi.value(2:Nx+2,2:Ny+1)-phi.value(1:Nx+1,2:Ny+1))./dx;
yvalue = (phi.value(2:Nx+1,2:Ny+2)-phi.value(2:Nx+1,1:Ny+1))./dy;
zvalue=[];
faceGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
