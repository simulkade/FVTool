function phiBC = cellBoundaryCylindrical3D(phi, BC)
% function phiBC = cellBoundary2D(MeshStructure, BC, phi)
% It creates the matrix of coefficient based on the BC structure provided
% by the user. It also generates the right hand side vector of the linear
% system of equations
%
% SYNOPSIS:
%   phiBC = cellBoundaryCylindrical3D(phi, BC)
%
% PARAMETERS:
%   MeshStructure: a mesh structure created by buildMesh* functions
%   BC: boundary condition structure created by createBC function
%   phi: cell variable created by createCellVariable
%
% RETURNS:
%   phiBC: a cell variable including the values of the ghost cells
%
% EXAMPLE:
%   m = buildMesh2D(3,4,1,1);
%   phi = createCellVariable(m,1);
%   bc = createBC(m);
%   phi_with_ghost = cellBoundary(phi, bc)
%
% SEE ALSO:
%


% extract data from the mesh structure
Nxyz = BC.domain.dims;
Nx = Nxyz(1); Ntetta = Nxyz(2); Nz = Nxyz(3);
dx_1 = BC.domain.cellsize.x(1);
dx_end = BC.domain.cellsize.x(end);
dtetta_1 = BC.domain.cellsize.y(1);
dtetta_end = BC.domain.cellsize.y(end);
dz_1 = BC.domain.cellsize.z(1);
dz_end = BC.domain.cellsize.z(end);
rp = repmat(BC.domain.cellcenters.x, 1, Nz);
% define the output matrix
phiBC = zeros(Nx+2, Ntetta+2, Nz+2);
phiBC(2:Nx+1, 2:Ntetta+1, 2:Nz+1) = phi;

% Assign values to the boundary values
if (BC.top.periodic==0) && (BC.bottom.periodic==0)
    % top boundary
    j=Ntetta+2;
    i = 2:Nx+1;
    k = 2:Nz+1;
    phiBC(i,j,k)= ...
        (BC.top.c-squeeze(phi(:,end,:)).*(-BC.top.a./(dtetta_end*rp)+BC.top.b/2))./(BC.top.a./(dtetta_end*rp)+BC.top.b/2);

    % Bottom boundary
    j=1;
    i = 2:Nx+1;
    k = 2:Nz+1;
    phiBC(i,j,k)= ...
        (BC.bottom.c-squeeze(phi(:,1,:)).*(BC.bottom.a./(dtetta_1*rp)+BC.bottom.b/2))./(-BC.bottom.a./(dtetta_1*rp)+BC.bottom.b/2);
else
    % top boundary
    j=Ntetta+2;
    i = 2:Nx+1;
    k = 2:Nz+1;
    phiBC(i,j,k)= phi(:,1,:);

    % Bottom boundary
    j=1;
    i = 2:Nx+1;
    k = 2:Nz+1;
    phiBC(i,j,k)= phi(:,end,:);
end

if (BC.left.periodic==0) && (BC.right.periodic==0)
    % Right boundary
    i = Nx+2;
    j = 2:Ntetta+1;
    k = 2:Nz+1;
    phiBC(i,j,k)= ...
        (BC.right.c-squeeze(phi(end,:,:)).*(-BC.right.a/dx_end+BC.right.b/2))./(BC.right.a/dx_end+BC.right.b/2);

    % Left boundary
    i = 1;
    j = 2:Ntetta+1;
    k = 2:Nz+1;
    phiBC(i,j,k)= ...
        (BC.left.c-squeeze(phi(1,:,:)).*(BC.left.a/dx_1+BC.left.b/2))./(-BC.left.a/dx_1+BC.left.b/2);
else
    % Right boundary
    i = Nx+2;
    j = 2:Ntetta+1;
    k = 2:Nz+1;
    phiBC(i,j,k)= phi(1,:,:);

    % Left boundary
    i = 1;
    j = 2:Ntetta+1;
    k = 2:Nz+1;
    phiBC(i,j,k)= phi(end,:,:);
end

if (BC.bottom.periodic==0) && (BC.top.periodic==0)
    % front boundary
    i = 2:Nx+1;
    j = 2:Ntetta+1;
    k = Nz+2;
    phiBC(i,j,k)= ...
        (BC.front.c-squeeze(phi(:,:,end)).*(-BC.front.a/dz_end+BC.front.b/2))./(BC.front.a/dz_end+BC.front.b/2);

    % back boundary
    i = 2:Nx+1;
    j = 2:Ntetta+1;
    k = 1;
    phiBC(i,j,k)= ...
        (BC.back.c-squeeze(phi(:,:,1)).*(BC.back.a/dz_1+BC.back.b/2))./(-BC.back.a/dz_1+BC.back.b/2);
else
    % front boundary
    i = 2:Nx+1;
    j = 2:Ntetta+1;
    k = Nz+2;
    phiBC(i,j,k)= phi(:,:,1);

    % back boundary
    i = 2:Nx+1;
    j = 2:Ntetta+1;
    k = 1;
    phiBC(i,j,k)= phi(:,:,end);
end
