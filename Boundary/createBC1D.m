function BC = createBC1D(meshvar)
% function BC = createBC1D(meshvar)
% Creates a boundary condition structure from a mesh structure
%
% SYNOPSIS:
%   BC = createBC1D(meshvar)
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

% Define the right, and left boundary conditions
% (default = Neumann, i.e., a = 1, b = 0, c = 0)

right.a = 1;
right.b = 0;
right.c = 0;
right.periodic = 0;

left.a = 1;
left.b = 0;
left.c = 0;
left.periodic = 0;

bottom=[];
top=[];
back=[];
front=[];

BC= BoundaryCondition(meshvar, left, right, bottom, top, back, front);
