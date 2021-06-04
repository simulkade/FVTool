function phiBC = cellBoundary(phi, BC)
% This function calculates the value of the boundary cells and add them
% to the variable phi. The function is used to add ghost cells to a matrix
% that specifies the cell values over a domain
%
% SYNOPSIS:
%   phiBC = cellBoundary(phi, BC)
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

% extract data from the mesh structure
d = BC.domain.dimension;

if (d ==1) || (d==1.5) || (d==1.8)
	phiBC = cellBoundary1D(phi, BC);
elseif (d == 2) || (d == 2.5)
	phiBC = cellBoundary2D(phi, BC);
elseif (d==2.8)
    phiBC = cellBoundaryRadial2D(phi, BC);
elseif (d==3)
    phiBC = cellBoundary3D(phi, BC);
elseif (d==3.2)
    phiBC = cellBoundaryCylindrical3D(phi, BC);
end
