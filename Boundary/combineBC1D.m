function [Mout, RHSout] = combineBC1D(BC, Meq, RHSeq)
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
Nx = BC.domain.dims(1);
dx_1 = BC.domain.cellsize.x(1);
dx_end= BC.domain.cellsize.x(end);
G = 1:Nx+2;

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
M(right) = M(right)-((BC.right.b/2 - BC.right.a/dx_end)./(BC.right.b/2 + BC.right.a/dx_end)).*M(rightE);
RHS(G(i-1)) = RHS(G(i-1))-M(rightE).*BC.right.c./(BC.right.b/2 + BC.right.a/dx_end);

% Left boundary
i = 1;
left = sub2ind(ms, G(i+1), G(i+1)); % left boundary cells
leftW = sub2ind(ms, G(i+1), G(i)); % west cells to left boundary cells
M(left) = M(left)-((BC.left.b/2 + BC.left.a/dx_1)./(BC.left.b/2 - BC.left.a/dx_1)).*M(leftW);
RHS(G(i+1)) = RHS(G(i+1))-M(leftW).*BC.left.c./(BC.left.b/2 - BC.left.a/dx_1);

Mout = M(G(2:end-1), G(2:end-1));
RHSout = RHS(reshape(G(2:end-1),Nx,1));
