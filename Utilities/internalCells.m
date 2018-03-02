function cellvar = internalCells(phi)
% this function extracts the internal cells from variable phi, which
% icludes ghost cells.
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

% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file

% check the size of the variable and the mesh dimension
d = phi.domain.dimension;
N = phi.domain.dims;

if (d ==1) || (d==1.5) || (d==1.8)
	cellvar = phi.value(2:N(1)+1);
elseif (d == 2) || (d == 2.5) || (d==2.8)
	cellvar = phi.value(2:N(1)+1, 2:N(2)+1);
elseif (d == 3) || (d==3.2)
    cellvar = phi.value(2:N(1)+1, 2:N(2)+1, 2:N(3)+1);
end
