function [Mout, RHSout] = combineBC1D(MeshStructure, BC, Meq, RHSeq)
%COMBINEBC This function combines the boundary condition equations with the
%main physical model equations, and delivers the matrix of coefficient and
%RHS to be solved for the internal cells.
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
dx = MeshStructure.cellsize;

% define the RHS column vector
ms = size(Meq);
M = Meq;
RHS = RHSeq;

% Assign values to the boundary condition matrix and the RHS vector based
% on the BC structure

% Right boundary
i=Nx+2;
right = sub2ind(ms, G(i-1), G(i-1)); % right boundary cells
rightE = sub2ind(ms, G(i-1), G(i)); % east cells to right boundary cells
M(right) = M(right)-((BC.right.b/2 - BC.right.a/dx)./(BC.right.b/2 + BC.right.a/dx)).*M(rightE);
RHS(G(i-1)) = RHS(G(i-1))-M(rightE).*BC.right.c./(BC.right.b/2 + BC.right.a/dx);

% Left boundary
i = 1;
left = sub2ind(ms, G(i+1), G(i+1)); % left boundary cells
leftW = sub2ind(ms, G(i+1), G(i)); % west cells to left boundary cells
M(left) = M(left)-((BC.left.b/2 + BC.left.a/dx)./(BC.left.b/2 - BC.left.a/dx)).*M(leftW);
RHS(G(i+1)) = RHS(G(i+1))-M(leftW).*BC.left.c./(BC.left.b/2 - BC.left.a/dx);

Mout = M(G(2:end-1), G(2:end-1));
RHSout = RHS(reshape(G(2:end-1),Nx,1));
