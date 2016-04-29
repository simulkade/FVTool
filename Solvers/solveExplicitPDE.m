function phi = solveExplicitPDE(MeshStructure, BC, RHS, dt, phi_old, varargin)
% SolveExplicitPDE solves for the new value of variable \phi in an explicit
% discretization scheme: \phi_new = \phi_old - dt/\alpha * RHS
% The code calculates the new values for the internal cells. Then it uses
% the boundary condition to calculate the values for the ghost cells

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

n = MeshStructure.dimension;
N = MeshStructure.numberofcells;

if (n>=2)
    phi = reshape(x, N+2);
else
    phi = reshape(x, [N+2 1]);
end


end
