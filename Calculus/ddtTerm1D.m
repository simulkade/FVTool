function ddt = ddtTerm1D(MeshStructure, dt, phi)
% function ddt = ddtTerm1D(MeshStructure, h, dt, phi)
% it returns the derivative of phi with respect to time in the vector ddt
% phi is a structure with elements phi.value and phi.Old
% mesh structure is used to map the (phi.value-phi.Old) matrix in the
% ddt vector
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

row_index = reshape(G(2:Nx+1),Nx,1); % main diagonal (only internal cells)

% define the RHS Vector
ddt = zeros(Nx+2,1);

% assign the values of the RHS vector
ddt(row_index) = reshape((phi.value(2:Nx+1)-phi.Old(2:Nx+1))/dt,Nx,1);
