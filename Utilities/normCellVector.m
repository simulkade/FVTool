function cellnorm = normCellVector(cellvec)
% this function calculates the second norm of a cell vector field
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

% check the size of the variable and the mesh dimension
d = cellvec.domain.dimension;

if (d ==1) || (d==1.5)
	cellnormval = abs(cellvec.xvalue);
elseif (d == 2) || (d == 2.5) || (d==2.8)
	cellnormval = realsqrt(cellvec.xvalue.*cellvec.xvalue+ ...
        cellvec.yvalue.*cellvec.yvalue);
elseif (d == 3) || (d==3.2)
    cellnormval = realsqrt(cellvec.xvalue.*cellvec.xvalue+ ...
        cellvec.yvalue.*cellvec.yvalue+ ...
        cellvec.zvalue.*cellvec.zvalue);
end
BC = createBC(cellvec.domain);
c=cellBoundary(cellnormval, BC);
cellnorm=CellVariable(cellvec.domain, c);
