function [BCMatrix, BCRHS] = boundaryCondition2D(BC)
% It creates the matrix of coefficient based on the BC structureprovided
% by the user. It also generates the right hand side vector of the linear
% system of equations
%
% SYNOPSIS:
%   [BCMatrix, BCRHS] = boundaryCondition2D(BC)
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

% extract data from the mesh structure
Nxy = BC.domain.dims;
Nx = Nxy(1); Ny = Nxy(2);
G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
dx_1 = BC.domain.cellsize.x(1);
dx_end = BC.domain.cellsize.x(end);
dy_1 = BC.domain.cellsize.y(1);
dy_end = BC.domain.cellsize.y(end);

% number of boundary nodes:
nb = 8*(Nx+Ny+2);

% define the vectors to be used for the creation of the sparse matrix
ii = zeros(nb,1);
jj = zeros(nb,1);
s = zeros(nb,1);

% define the RHS column vector
BCRHS = zeros((Nx+2)*(Ny+2), 1);

% assign value to the corner nodes (useless cells)
q = 1:4;
ii(q) = BC.domain.corners; jj(q) = BC.domain.corners;
s(q) = max(BC.top.b/2 + BC.top.a/dy_end); BCRHS(BC.domain.corners) = 0;

% Assign values to the boundary condition matrix and the RHS vector based
% on the BC structure

if (BC.top.periodic ==0) && (BC.bottom.periodic ==0)
    % top boundary
    j=Ny+2;
    i=2:Nx+1;
    q = q(end)+(1:Nx);
    ii(q) = G(i,j);  jj(q) = G(i,j);  s(q) = BC.top.b/2 + BC.top.a/dy_end;
    q = q(end)+(1:Nx);
    ii(q) = G(i,j);  jj(q) = G(i,j-1); s(q) = BC.top.b/2 - BC.top.a/dy_end;
    BCRHS(G(i,j)) = BC.top.c;

    % Bottom boundary
    j=1;
    i=2:Nx+1;
    q = q(end)+(1:Nx);
    ii(q) = G(i,j);  jj(q) = G(i,j+1);  s(q) = -(BC.bottom.b/2 + BC.bottom.a/dy_1);
    q = q(end)+(1:Nx);
    ii(q) = G(i,j);  jj(q) = G(i,j); s(q) = -(BC.bottom.b/2 - BC.bottom.a/dy_1);
    BCRHS(G(i,j)) = -(BC.bottom.c);
elseif (BC.top.periodic ==1) || (BC.bottom.periodic ==1) % periodic boundary
    % top boundary
    j=Ny+2;
    i=2:Nx+1;
    q = q(end)+(1:Nx);
    ii(q) = G(i,j);  jj(q) = G(i,j);  s(q) = 1;
    q = q(end)+(1:Nx);
    ii(q) = G(i,j);  jj(q) = G(i,j-1);  s(q) = -1;
    q = q(end)+(1:Nx);
    ii(q) = G(i,j);  jj(q) = G(i,1); s(q) = dy_end/dy_1;
    q = q(end)+(1:Nx);
    ii(q) = G(i,j);  jj(q) = G(i,2); s(q) = -dy_end/dy_1;
    BCRHS(G(i,j)) = 0;

    % Bottom boundary
    j=1;
    i=2:Nx+1;
    q = q(end)+(1:Nx);
    ii(q) = G(i,j);  jj(q) = G(i,j);  s(q) = 1;
    q = q(end)+(1:Nx);
    ii(q) = G(i,j);  jj(q) = G(i,j+1);  s(q) = 1;
    q = q(end)+(1:Nx);
    ii(q) = G(i,j);  jj(q) = G(i,Ny+1); s(q) = -1;
    q = q(end)+(1:Nx);
    ii(q) = G(i,j);  jj(q) = G(i,Ny+2); s(q) = -1;
    BCRHS(G(i,j)) = 0;
end


if (BC.right.periodic == 0) && (BC.left.periodic == 0)
    % Right boundary
    i=Nx+2;
    j=2:Ny+1;
    q = q(end)+(1:Ny);
    ii(q) = G(i,j);  jj(q) = G(i,j);  s(q) = BC.right.b/2 + BC.right.a/dx_end;
    q = q(end)+(1:Ny);
    ii(q) = G(i,j);  jj(q) = G(i-1,j); s(q) = BC.right.b/2 - BC.right.a/dx_end;
    BCRHS(G(i,j)) = BC.right.c;

    % Left boundary
    i = 1;
    j=2:Ny+1;
    q = q(end)+(1:Ny);
    ii(q) = G(i,j);  jj(q) = G(i+1,j);  s(q) = -(BC.left.b/2 + BC.left.a/dx_1);
    q = q(end)+(1:Ny);
    ii(q) = G(i,j);  jj(q) = G(i,j); s(q) = -(BC.left.b/2 - BC.left.a/dx_1);
    BCRHS(G(i,j)) = -(BC.left.c);
elseif (BC.right.periodic == 1) || (BC.left.periodic == 1) % periodic boundary
    % Right boundary
    i=Nx+2;
    j=2:Ny+1;
    q = q(end)+(1:Ny);
    ii(q) = G(i,j);  jj(q) = G(i,j);  s(q) = 1;
    q = q(end)+(1:Ny);
    ii(q) = G(i,j);  jj(q) = G(i-1,j);  s(q) = -1;
    q = q(end)+(1:Ny);
    ii(q) = G(i,j);  jj(q) = G(1,j); s(q) = dx_end/dx_1;
    q = q(end)+(1:Ny);
    ii(q) = G(i,j);  jj(q) = G(2,j); s(q) = -dx_end/dx_1;
    BCRHS(G(i,j)) = 0;

    % Left boundary
    i = 1;
    j=2:Ny+1;
    q = q(end)+(1:Ny);
    ii(q) = G(i,j);  jj(q) = G(i,j);  s(q) = 1;
    q = q(end)+(1:Ny);
    ii(q) = G(i,j);  jj(q) = G(i+1,j);  s(q) = 1;
    q = q(end)+(1:Ny);
    ii(q) = G(i,j);  jj(q) = G(Nx+1,j); s(q) = -1;
    q = q(end)+(1:Ny);
    ii(q) = G(i,j);  jj(q) = G(Nx+2,j); s(q) = -1;
    BCRHS(G(i,j)) = 0;
end
% Build the sparse matrix of the boundary conditions
BCMatrix = sparse(ii(1:q(end)), jj(1:q(end)), s(1:q(end)), ...
    (Nx+2)*(Ny+2), (Nx+2)*(Ny+2));
