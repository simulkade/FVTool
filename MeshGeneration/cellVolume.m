function cellvol = cellVolume(meshvar)
% cellvol = cellVolume(meshvar)
% returns the volume of each cell as a cell variable
% SYNOPSIS:
%   cellvol = cellVolume(meshvar)
%
% PARAMETERS:
%   MeshStructure: a mesh structure created by buildMesh* functions
%
% RETURNS:
%   cellvol: a (1D, 2D, or 3D) matrix depending on the mesh size
%
% EXAMPLE:
%   m = createMesh2D(3,4, 1.0, 2.0); % creates a mesh
%   cell_vol=cellVolume(m);
%
% SEE ALSO:
%     createFaceVariable, createBC, buildMesh1D,
%     buildMesh2D, buildMesh3D,
%     buildMeshCylindrical1D, buildMeshCylindrical2D,
%     cellBoundary, combineBC

% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file

% check the size of the variable and the mesh dimension
dim = meshvar.dimension;
BC = createBC(meshvar);
switch dim
    case 1
        c=meshvar.cellsize.x(2:end-1);
    case 1.5
        c=2.0*pi()*meshvar.cellsize.x(2:end-1).*meshvar.cellcenters.x;
    case 2
        c=meshvar.cellsize.x(2:end-1)*meshvar.cellsize.y(2:end-1)';
    case 2.5 % cylindrical
        c=2.0*pi()*meshvar.cellcenters.x.*meshvar.cellsize.x(2:end-1)*meshvar.cellsize.y(2:end-1)';
    case 2.8 % radial
        c=meshvar.cellcenters.x.*meshvar.cellsize.x(2:end-1)*meshvar.cellsize.y(2:end-1)';
    case 3
        Nx = meshvar.dims(1);
        Ny = meshvar.dims(2);
        Nz = meshvar.dims(3);
        DXp = repmat(meshvar.cellsize.x(2:end-1), 1, Ny, Nz);
        DYp = repmat(meshvar.cellsize.y(2:end-1)', Nx, 1, Nz);
        DZ = zeros(1,1,Nz+2);
        DZ(1,1,:) = meshvar.cellsize.z;
        DZp=repmat(DZ(1,1,2:end-1), Nx, Ny, 1);
        c=DXp.*DYp.*DZp;
    case 3.2
        N = meshvar.dims;
        Nr = N(1); Ntetta=N(2); Nz = N(3);
        rp = repmat(meshvar.cellcenters.x, 1, Ntetta, Nz);
        DRp = repmat(meshvar.cellsize.x(2:end-1), 1, Ntetta, Nz);
        DTHETAp = repmat(meshvar.cellsize.y(2:end-1)', Nr, 1, Nz);
        DZ = zeros(1,1,Nz+2);
        DZ(1,1,:) = meshvar.cellsize.z;
        DZp=repmat(DZ(1,1,2:end-1), Nr, Ntetta, 1);
        c=rp.*DRp.*DTHETAp.*DZp;
end
cellvol= createCellVariable(meshvar, c, BC);
