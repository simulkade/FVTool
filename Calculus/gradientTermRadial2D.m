function faceGrad = gradientTermRadial2D(phi)
% this function calculates the gradient of a variable in x and y direction
% it checks for the availability of the ghost variables and use them, otherwise
% estimate them, assuming a zero gradient on the boundaries
%
% SYNOPSIS:
%   faceGrad = gradientTermRadial2D(phi)
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
Nr = phi.domain.dims(1);
Ntheta = phi.domain.dims(2);
DR = repmat(phi.domain.cellsize.x, 1, Ntheta);
DTHETA = repmat(phi.domain.cellsize.y', Nr, 1);
dr = 0.5*(DR(1:end-1,:)+DR(2:end,:));
dtetta = 0.5*(DTHETA(:,1:end-1)+DTHETA(:,2:end));
rp = repmat(phi.domain.cellcenters.x, 1, Ntheta+1);

% in this case, ghost cells have values
xvalue = (phi.value(2:Nr+2,2:Ntheta+1)-phi.value(1:Nr+1,2:Ntheta+1))./dr;
yvalue = (phi.value(2:Nr+1,2:Ntheta+2)-phi.value(2:Nr+1,1:Ntheta+1))./(dtetta.*rp);
zvalue=[];
faceGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
