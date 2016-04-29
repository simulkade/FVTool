function RHS = constantSourceTerm2D(phi)
% RHS vector for a Explicit source term
%
% SYNOPSIS:
%   RHS = constantSourceTerm2D(phi)
%
% PARAMETERS:
%   phi: Cell Variable
%
% RETURNS:
%   RHS: vector
%
% EXAMPLE:
%
% SEE ALSO:
%

% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file

% extract data from the mesh structure
Nxy = phi.domain.dims;
Nx = Nxy(1); Ny = Nxy(2);
G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);

% rearrange the matrix of k and build the sparse matrix for internal cells
row_index = reshape(G(2:Nx+1,2:Ny+1),Nx*Ny,1); % main diagonal (only internal cells)

% define the RHS Vector
RHS = zeros((Nx+2)*(Ny+2),1);

% assign the values of the RHS vector
RHS(row_index) = reshape(phi.value(2:Nx+1,2:Ny+1),Nx*Ny,1);
