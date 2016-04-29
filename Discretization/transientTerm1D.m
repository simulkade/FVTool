function [M, RHS] = transientTerm1D(MeshStructure, h, dt, phi)
% Matrix of coefficients and the RHS vector for a transient term
% h \partial_t \phi
% h is a matrix of size [m,1]
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
Nx = MeshStructure.numberofcells;

% rearrange the matrix of k and build the sparse matrix for internal cells
row_index = reshape(G(2:Nx+1),Nx,1); % main diagonal (only internal cells)
AP_diag = reshape(h/dt,Nx,1);
M = sparse(row_index, row_index, AP_diag, Nx+2, Nx+2);

% extract the previous value of \phi (including the ghost cells)
phi_old = phi.Old;

% define the RHS Vector
RHS = zeros(Nx+2,1);

% assign the values of the RHS vector
RHS(row_index) = reshape(h.*phi_old(2:Nx+1)/dt,Nx,1);
