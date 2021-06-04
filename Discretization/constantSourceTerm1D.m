function RHS = constantSourceTerm1D(phi)
% RHS vector for a Explicit source term
%
% SYNOPSIS:
%   RHS = constantSourceTerm1D(phi)
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
Nx = phi.domain.dims(1);
G = 1:Nx+2;

% rearrange the matrix of k and build the sparse matrix for internal cells
row_index = reshape(G(2:Nx+1),Nx,1); % main diagonal (only internal cells)

% define the RHS Vector
RHS = zeros(Nx+2,1);

% assign the values of the RHS vector
RHS(row_index) = reshape(phi.value(2:Nx+1),Nx,1);
