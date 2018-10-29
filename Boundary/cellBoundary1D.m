function phiBC = cellBoundary1D(phi, BC)
% This function calculates the value of the boundary cells and add them
% to the variable phi of size (1 .. Nx)
% the output includes the boundary is of the size (1..Nx+2)
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

% extract data from the mesh structure
% Nx = MeshStructure.numberofcells;
dx_1 = BC.domain.cellsize.x(1);
dx_end= BC.domain.cellsize.x(end);

% boundary condition (a d\phi/dx + b \phi = c, a column vector of [d a])
% a (phi(i)-phi(i-1))/dx + b (phi(i)+phi(i-1))/2 = c
% phi(i) (a/dx+b/2) + phi(i-1) (-a/dx+b/2) = c
% Right boundary, i=m+2
% phi(i) (a/dx+b/2) = c- phi(i-1) (-a/dx+b/2)
% Left boundary, i=2
%  phi(i-1) (-a/dx+b/2) = c - phi(i) (a/dx+b/2)
% define the new phi
if (BC.left.periodic==0) && (BC.right.periodic==0)
    phiBC = [(BC.left.c-phi(1)* ...
        (BC.left.a/dx_1+BC.left.b/2))/(-BC.left.a/dx_1+BC.left.b/2); phi; ...
        (BC.right.c-phi(end)* ...
        (-BC.right.a/dx_end+BC.right.b/2))/(BC.right.a/dx_end+BC.right.b/2)];
else
    phiBC = [phi(end); phi; phi(1)];
end
