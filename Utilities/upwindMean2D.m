function phiFaceAverage = upwindMean2D(phi, u)
% This function gets the value of the field variable phi defined
% over the MeshStructure and calculates the upwind average on
% the cell faces, based on the direction of the velocity vector for a uniform mesh.
%
% SYNOPSIS:
%   phiFaceAverage = upwindMean2D(phi, u)
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

% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue;
uy = u.yvalue;

% check the size of the variable and the mesh dimension
Nxy = phi.domain.dims;
Nx = Nxy(1); Ny = Nxy(2);

% assign to a temp variable for boundary corrections
phi_tmp = phi.value;

% correct the value of phi at the boundary (calculation trick)
% assign the value of the left boundary to the left ghost cells
phi_tmp(1,:) = (phi.value(1,:)+phi.value(2,:))/2;
% assign the value of the right boundary to the right ghost cells
phi_tmp(end,:) = (phi.value(end,:)+phi.value(end-1,:))/2;
% assign the value of the bottom boundary to the bottom ghost cells
phi_tmp(:,1) = (phi.value(:,1)+phi.value(:,2))/2;
% assign the value of the top boundary to the top ghost cells
phi_tmp(:,end) = (phi.value(:,end)+phi.value(:,end-1))/2;

% calculate the average value
xvalue = (ux>0).*phi_tmp(1:Nx+1,2:Ny+1)+ ...
                        (ux<0).*phi_tmp(2:Nx+2,2:Ny+1)+ ...
                        0.5*(ux==0).*(phi.value(1:Nx+1,2:Ny+1)+phi.value(2:Nx+2,2:Ny+1));
yvalue = (uy>0).*phi_tmp(2:Nx+1,1:Ny+1)+ ...
                        (uy<0).*phi_tmp(2:Nx+1,2:Ny+2)+ ...
                        0.5*(uy==0).*(phi.value(2:Nx+1,1:Ny+1)+phi.value(2:Nx+1,2:Ny+2));
phiFaceAverage=FaceVariable(phi.domain, xvalue, yvalue, []);
