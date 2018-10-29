function BC = createBC1D(meshvar)
% function BC = createBC1D(meshvar)
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
