function cellvar = createCellVariable(meshvar, cellval, varargin)
% cellvar = createCellVariable(MeshStructure, cellval)
% cellvar = createCellVariable(MeshStructure, cellval, BC)
% this function creates a cell variable based on the dimension of the
% problem. The cell variable does not include the ghost cells, unless
% the boundary condition is supplied as an input. The value of
% the ghost cells can be added to the cell variable using the
% cellBoundary function or by supplying the BC structure as the third
% argument of this function.
% SYNOPSIS:
%   cellvar = createCellVariable(MeshStructure, cellval)
%   cellvar = createCellVariable(MeshStructure, cellval, BC)
%
% PARAMETERS:
%   MeshStructure: a mesh structure created by buildMesh* functions
%   cellValue: the value that will be assigned to the internal cells
%   BC (optional): boundary condition structure created by createBC
%   function
%
% RETURNS:
%   cellvar: a (1D, 2D, or 3D) matrix depending on the mesh size
%
% EXAMPLE:
%   m = createMesh2D(3,4, 1.0, 2.0); % creates a mesh
%   bc = createBC(m); % creates a Neumann boundary condition
%   phi_val = 2.0; % value to be assigned to all cells
%   phi = createCellVariable(m, phi_val); % value assigned to the internal cells with Neumann BC
%   phi_with_ghost = createCellVariable(m, phi_val, bc); % value assigned to the internal and ghost cells
%
% SEE ALSO:
%     createFaceVariable, createBC, createMesh1D,
%     createMesh2D, createMesh3D, createMeshRadial2D
%     createMeshCylindrical1D, createMeshCylindrical2D,
%     createMeshCylindrical3D, cellBoundary, combineBC

% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file

% check the size of the variable and the mesh dimension
dim = meshvar.dims;
if nargin<3
    BC = createBC(meshvar);
else
    BC=varargin{1};
end
if numel(cellval)==1
    c=cellval*ones(dim);
elseif prod(size(cellval)==dim)
    c=cellval;
else
    c= zeros(dim);
end
c=cellBoundary(c, BC);
cellvar= CellVariable(meshvar, c);
