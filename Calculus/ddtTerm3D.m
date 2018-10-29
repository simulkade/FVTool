function ddt = ddtTerm3D(MeshStructure, dt, phi)
% Matrix of coefficients and the RHS vector for a transient term
% h \partial_t \phi
% h is a matrix of size [m, n]
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

% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file

% extract data from the mesh structure
G = MeshStructure.numbering;
Nxyz = MeshStructure.numberofcells;
Nx = Nxyz(1); Ny = Nxyz(2); Nz = Nxyz(3);

row_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),Nx*Ny*Nz,1); % main diagonal (only internal cells)

% define the RHS Vector
ddt = zeros((Nx+2)*(Ny+2)*(Nz+2),1);

% assign the values of the RHS vector
ddt(row_index) = ...
    reshape((phi.value(2:Nx+1,2:Ny+1,2:Nz+1)-phi.Old(2:Nx+1,2:Ny+1,2:Nz+1))/dt,Nx*Ny*Nz,1);
