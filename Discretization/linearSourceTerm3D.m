function M = linearSourceTerm3D(k)
% Matrix of coefficients for a linear source term in the form of k \phi
%
%
% SYNOPSIS:
%   M = linearSourceTerm3D(k)
%
% PARAMETERS:
%   k: CellVariable
%
% RETURNS:
%   M: sparse matrix
%
% EXAMPLE:
%
% SEE ALSO:
%

% extract data from the mesh structure
Nxyz = k.domain.dims;
Nx = Nxyz(1); Ny = Nxyz(2); Nz = Nxyz(3);
G=reshape(1:(Nx+2)*(Ny+2)*(Nz+2), Nx+2, Ny+2, Nz+2);

% rearrange the matrix of k and build the sparse matrix for internal cells
row_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),Nx*Ny*Nz,1); % main diagonal (only internal cells)
AP_diag = reshape(k.value(2:Nx+1, 2:Ny+1, 2:Nz+1),Nx*Ny*Nz,1);
M = sparse(row_index, row_index, AP_diag, ...
    (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2));
