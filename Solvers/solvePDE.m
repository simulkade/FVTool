function phi = solvePDE(MS, M, RHS, varargin)
%SOLVEPDE solves the linear system M x \phi = RHS and returns the value of
%\phi reshaped based on the structure of the MeshStructure variable. The
%default solver is the matlab '\' linear solver
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

x = M\RHS;

n = MS.dimension;
N = MS.dims;

if (n>=2)
    phival = reshape(x, N+2);
else
    phival = reshape(x, [N(1)+2 1]);
end

phi=CellVariable(MS, phival);

end
