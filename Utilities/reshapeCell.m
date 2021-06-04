function cellvar = reshapeCell(MS, phi)
% this function reshapes a vetorized cell variable to its domain shape
% matrix based on the mesh structure data; it is assumed that the phi
% includes the ghost cell data as well.
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
%     cellBoundary

% Written by Ali A. Eftekhari
% See the license file

% check the size of the variable and the mesh dimension
d = MS.dimension;
N = MS.dims;

if (d ==1) || (d==1.5) || (d==1.8)
	cellvar = reshape(phi, N(1)+2, 1);
elseif (d == 2) || (d == 2.5) || (d==2.8)
	cellvar = reshape(phi, N(1)+2, N(2)+2);
elseif (d == 3) || (d==3.2)
    cellvar = reshape(phi, N(1)+2, N(2)+2, N(3)+2);
end
