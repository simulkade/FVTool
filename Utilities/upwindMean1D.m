function phiFaceAverage = upwindMean1D(phi, u)
% This function gets the value of the field variable phi defined
% over the MeshStructure and calculates the upwind average on
% the cell faces, based on the direction of the velocity vector for a uniform mesh.
%
% SYNOPSIS:
%   phiFaceAverage = upwindMean1D(phi, u)
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

% Written by Ali A. Eftekhari
% See the license file

% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue;

% check the size of the variable and the mesh dimension
Nx = phi.domain.dims(1);

% assign to a temp variable for boundary corrections
phi_tmp = phi.value;

% correct the value of phi at the boundary (calculation trick)
% assign the value of the left boundary to the left ghost cell
phi_tmp(1) = (phi.value(1)+phi.value(2))/2;
% assign the value of the right boundary to the right ghost cell
phi_tmp(end) = (phi.value(end)+phi.value(end-1))/2;

% calculate the average value
xvalue = (ux>0).*phi_tmp(1:Nx+1)+ ...
                        (ux<0).*phi_tmp(2:Nx+2)+ ...
                        0.5*(ux==0).*(phi.value(1:Nx+1)+phi.value(2:Nx+2));
phiFaceAverage=FaceVariable(phi.domain, xvalue, [], []);
