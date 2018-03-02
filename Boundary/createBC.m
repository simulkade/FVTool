function BC = createBC(meshvar)
% function BC = createBC(meshvar)
% Creates a boundary condition structure from a mesh structure
% for a 2D mesh. The boundary conditions on all boundaries are Neumann;
% The index of each boundary condition is defined as:
%	1:	Dirichlet
%	2:	Neumann
%	3:	Mixed
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

d = meshvar.dimension;
if (d == 1) || (d==1.5) || (d==1.8)
	BC = createBC1D(meshvar);
elseif (d == 2) || (d == 2.5) || (d == 2.8)
	BC = createBC2D(meshvar);
elseif (d == 3) || (d==3.2)
    BC = createBC3D(meshvar);
end
