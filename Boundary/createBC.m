function BC = createBC(meshvar)
% function BC = createBC(meshvar)
% Creates a boundary condition structure from a mesh structure
% The boundary conditions on all boundaries are Neumann by default
% and can be altered later bu the user; see examples
%
% SYNOPSIS:
%	BC = createBC(meshvar)
%
% PARAMETERS:
%
%
% RETURNS:
%
%
% EXAMPLE:
%	L = 1.0; % length of a 1D domain
%   Nx = 5; % Number of grids in x direction
%   m = createMesh1D(Nx, L);
%   BC = createBC(m); % all Neumann boundaries
%
% SEE ALSO:
%

d = meshvar.dimension;
if (d == 1) || (d==1.5) || (d==1.8)
	BC = createBC1D(meshvar);
elseif (d == 2) || (d == 2.5) || (d == 2.8)
	BC = createBC2D(meshvar);
elseif (d == 3) || (d==3.2)
    BC = createBC3D(meshvar);
end
