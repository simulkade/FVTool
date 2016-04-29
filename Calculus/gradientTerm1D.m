function faceGrad = gradientTerm1D(phi)
% this function calculates the gradient of a variable in x direction
% it checks for the availability of the ghost variables and use them, otherwise
% estimate them, assuming a zero gradient on the boundaries
%
% SYNOPSIS:
%   faceGrad = gradientTerm1D(phi)
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
DX = phi.domain.cellsize.x;
dx = 0.5*(DX(1:end-1)+DX(2:end));

% in this case, ghost cells have values
xvalue = (phi.value(2:Nx+2)-phi.value(1:Nx+1))./dx;
yvalue=[];
zvalue=[];

faceGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
