function [BCMatrix, BCRHS] = boundaryCondition(BC)
% creates the matrix of coefficient and RHS vector
%
% SYNOPSIS:
%   [BCMatrix, BCRHS] = boundaryCondition(MeshStructure, BC)
%
% PARAMETERS:
%   MeshStructure  - a mesh structure created by buildMesh* functions
%   BC             - Right hand side values of the boundary condition
%   equations
%
% RETURNS:
%   BCMatrix  - an square sparse matrix
%   BCRHS     - a column vector of values
%
% EXAMPLE:
%   L = 1.0; % length of a 1D domain
%   Nx = 5; % Number of grids in x direction
%   m = createMesh1D(Nx, L);
%   BC = createBC(m); % all Neumann boundaries
%   [Mbc, RHSbc] = boundaryCondition(BC);
%   spy(Mbc); % see inside the boundary matrix of coeffiecients
% SEE ALSO:
%     createBC, createMesh1D, createMesh2D, createMesh3D,
%     createMeshCylindrical1D, createMeshCylindrical2D,
%     createMeshRadial2D, createMeshCylindrical3D,
%     cellBoundary, combineBC, createCellVariable

% Copyright (c) 2012-2019 Ali Akbar Eftekhari
% See the license file

d = BC.domain.dimension;
if (d ==1) || (d==1.5) || (d==1.8)
	[BCMatrix, BCRHS] = boundaryCondition1D(BC);
elseif (d == 2) || (d == 2.5)
	[BCMatrix, BCRHS] = boundaryCondition2D(BC);
elseif (d == 2.8)
    [BCMatrix, BCRHS] = boundaryConditionRadial2D(BC);
elseif (d == 3)
    [BCMatrix, BCRHS] = boundaryCondition3D(BC);
elseif (d == 3.2)
    [BCMatrix, BCRHS] = boundaryConditionCylindrical3D(BC);
end
