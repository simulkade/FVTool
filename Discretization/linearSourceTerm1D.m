function M = linearSourceTerm1D(k)
% Matrix of coefficients for a linear source term in the form of k \phi
% in 1D cartesian grid.
% k is a matrix of size [m,1], i.e., internal cells
%
% SYNOPSIS:
%
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
Nx = k.domain.dims(1);
G = 1:Nx+2;

% rearrange the matrix of k and build the sparse matrix for internal cells
row_index = reshape(G(2:Nx+1),Nx,1); % main diagonal (only internal cells)
AP_diag = reshape(k.value(2:Nx+1),Nx,1);
M = sparse(row_index, row_index, AP_diag, Nx+2, Nx+2);
