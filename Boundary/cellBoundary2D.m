function phiBC = cellBoundary2D(phi, BC)
% function phiBC = cellBoundary2D(BC, phi)
% It creates the matrix of coefficient based on the BC structure provided
% by the user. It also generates the right hand side vector of the linear
% system of equations
%
% SYNOPSIS:
%   phiBC = cellBoundary2D(phi, BC)
%
% PARAMETERS:
%   BC: boundary condition structure created by createBC function
%   phi: cell variable created by createCellVariable
%
% RETURNS:
%   phiBC: a cell variable including the values of the ghost cells
%
% EXAMPLE:
%   m = createMesh2D(3,4,1,1);
%   phi = createCellVariable(m,1);
%   bc = createBC(m);
%   phi_with_ghost = cellBoundary(phi,bc)
%
% SEE ALSO:
%

% Note: I use a for loop here for more readability of the code!

% extract data from the mesh structure
Nxy = BC.domain.dims;
Nx = Nxy(1); Ny = Nxy(2);
dx_1 = BC.domain.cellsize.x(1);
dx_end = BC.domain.cellsize.x(end);
dy_1 = BC.domain.cellsize.y(1);
dy_end = BC.domain.cellsize.y(end);

% define the output matrix
phiBC = zeros(Nx+2, Ny+2);
phiBC(2:Nx+1, 2:Ny+1) = phi;

% Assign values to the boundary values
if (BC.top.periodic==0) && (BC.bottom.periodic==0)
    % top boundary
    j=Ny+2;
    i = 2:Nx+1;
    phiBC(i,j)= ...
        (BC.top.c-phi(:,end).*(-BC.top.a/dy_end+BC.top.b/2))./(BC.top.a/dy_end+BC.top.b/2);

    % Bottom boundary
    j=1;
    i = 2:Nx+1;
    phiBC(i,j)= ...
        (BC.bottom.c-phi(:,1).*(BC.bottom.a/dy_1+BC.bottom.b/2))./(-BC.bottom.a/dy_1+BC.bottom.b/2);
else
    % top boundary
    j=Ny+2;
    i = 2:Nx+1;
    phiBC(i,j)= phi(:,1);

    % Bottom boundary
    j=1;
    i = 2:Nx+1;
    phiBC(i,j)= phi(:,end);
end

if (BC.left.periodic==0) && (BC.right.periodic==0)
    % Right boundary
    i = Nx+2;
    j = 2:Ny+1;
    phiBC(i,j)= ...
        (BC.right.c-phi(end,:).*(-BC.right.a/dx_end+BC.right.b/2))./(BC.right.a/dx_end+BC.right.b/2);

    % Left boundary
    i = 1;
    j = 2:Ny+1;
    phiBC(i,j)= ...
        (BC.left.c-phi(1,:).*(BC.left.a/dx_1+BC.left.b/2))./(-BC.left.a/dx_1+BC.left.b/2);
else
    % Right boundary
    i = Nx+2;
    j = 2:Ny+1;
    phiBC(i,j)= phi(1,:);

    % Left boundary
    i = 1;
    j = 2:Ny+1;
    phiBC(i,j)= phi(end,:);
end
