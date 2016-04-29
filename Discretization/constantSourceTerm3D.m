function RHS = constantSourceTerm3D(phi)
% RHS vector for a Explicit source term
%
% SYNOPSIS:
%   RHS = constantSourceTerm3D(phi)
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

% extract data from the mesh structure
Nxyz = phi.domain.dims;
Nx = Nxyz(1); Ny = Nxyz(2); Nz = Nxyz(3);
G=reshape(1:(Nx+2)*(Ny+2)*(Nz+2), Nx+2, Ny+2, Nz+2);

% rearrange the matrix of k and build the sparse matrix for internal cells
row_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),Nx*Ny*Nz,1); % main diagonal (only internal cells)

% define the RHS Vector
RHS = zeros((Nx+2)*(Ny+2)*(Nz+2),1);

% assign the values of the RHS vector
RHS(row_index) = reshape(phi.value(2:Nx+1, 2:Ny+1, 2:Nz+1),Nx*Ny*Nz,1);
