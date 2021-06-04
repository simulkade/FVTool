function M = linearSourceTerm2D(k)
% Matrix of coefficients for a linear source term in the form of k \phi
%
%
% SYNOPSIS:
%   M = linearSourceTerm2D(k)
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
Nxy = k.domain.dims;
Nx = Nxy(1); Ny = Nxy(2);
G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);

% rearrange the matrix of k and build the sparse matrix for internal cells
row_index = reshape(G(2:Nx+1,2:Ny+1),Nx*Ny,1); % main diagonal (only internal cells)
AP_diag = reshape(k.value(2:Nx+1,2:Ny+1),Nx*Ny,1);
M = sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2), (Nx+2)*(Ny+2));
