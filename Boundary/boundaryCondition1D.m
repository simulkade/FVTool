function [BCMatrix, BCRHS] = boundaryCondition1D(BC)
% It creates the matrix of coefficient based on the BC structure provided
% by the user. It also generates the right hand side vector of the linear
% system of equations
%
% SYNOPSIS:
%   [BCMatrix, BCRHS] = boundaryCondition1D(BC)
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
Nx = BC.domain.dims(1);
dx_1 = BC.domain.cellsize.x(1);
dx_end= BC.domain.cellsize.x(end);
G = [1:Nx+2];

% number of boundary nodes:
nb = 8;

% define the vectors to be used for the creation of the sparse matrix
ii = zeros(nb,1);
jj = zeros(nb,1);
s = zeros(nb,1);

% define the RHS column vector
BCRHS = zeros(Nx+2, 1);
q = 0;
% Assign values to the boundary condition matrix and the RHS vector based
% on the BC structure
if (BC.right.periodic == 0) || (BC.left.periodic == 0) % non-periodic boundary condition
    % Right boundary
    i = Nx+2;
    q=q+1;
    ii(q) = G(i);  jj(q) = G(i);  s(q) = BC.right.b/2 + BC.right.a/dx_end;
    q=q+1;
    ii(q) = G(i);  jj(q) = G(i-1); s(q) = BC.right.b/2 - BC.right.a/dx_end;
    BCRHS(G(i)) = BC.right.c;

    % Left boundary
    i = 1;
    q=q+1;
    ii(q) = G(i);  jj(q) = G(i+1);  s(q) = -(BC.left.b/2 + BC.left.a/dx_1);
    q=q+1;
    ii(q) = G(i);  jj(q) = G(i); s(q) = -(BC.left.b/2 - BC.left.a/dx_1);
    BCRHS(G(i)) = -(BC.left.c);
elseif (BC.right.periodic == 1) || (BC.left.periodic == 1) % periodic boundary condition
    % Right boundary
    i = Nx+2;
    q=q+1;
    ii(q) = G(i);  jj(q) = G(i);  s(q) = 1;
    q=q+1;
    ii(q) = G(i);  jj(q) = G(i-1);  s(q) = -1;
    q=q+1;
    ii(q) = G(i);  jj(q) = G(1); s(q) = dx_end/dx_1;
    q=q+1;
    ii(q) = G(i);  jj(q) = G(2); s(q) = -dx_end/dx_1;
    BCRHS(G(i)) = 0;
    % Left boundary
    i = 1;
    q=q+1;
    ii(q) = G(i);  jj(q) = G(i);  s(q) = 1;
    q=q+1;
    ii(q) = G(i);  jj(q) = G(2);  s(q) = 1;
    q=q+1;
    ii(q) = G(i);  jj(q) = G(Nx+1); s(q) = -1;
    q=q+1;
    ii(q) = G(i);  jj(q) = G(Nx+2); s(q) = -1;
    BCRHS(G(i)) = 0;
end

% Build the sparse matrix of the boundary conditions
BCMatrix = sparse(ii(1:q), jj(1:q), s(1:q), Nx+2, Nx+2);
