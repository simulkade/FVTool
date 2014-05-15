function [Mout, RHSout] = combineBC2D(MeshStructure, BC, Meq, RHSeq)
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

%{
Copyright (c) 2012, 2013, Ali Akbar Eftekhari
All rights reserved.

Redistribution and use in source and binary forms, with or 
without modification, are permitted provided that the following 
conditions are met:

    *   Redistributions of source code must retain the above copyright notice, 
        this list of conditions and the following disclaimer.
    *   Redistributions in binary form must reproduce the above 
        copyright notice, this list of conditions and the following 
        disclaimer in the documentation and/or other materials provided 
        with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%}

% extract data from the mesh structure
G = MeshStructure.numbering;
Nxy = MeshStructure.numberofcells;
Nx = Nxy(1); Ny = Nxy(2);
dxdy = MeshStructure.cellsize;
dx = dxdy(1); dy = dxdy(2);

% define the RHS column vector
ms = size(Meq);
M = Meq;
RHS = RHSeq;

% Assign values to the boundary condition matrix and the RHS vector based
% on the BC structure
% top boundary
j=Ny+2;
i=2:Nx+1;
top = reshape(sub2ind(ms, G(i,j-1), G(i,j-1)), Nx,1); % top boundary cells
topN = reshape(sub2ind(ms, G(i,j-1), G(i,j)), Nx, 1); % north cells to top boundary cells
M(top) = M(top)-((BC.top.b/2 - BC.top.a/dy)./(BC.top.b/2 + BC.top.a/dy)).*M(topN);
RHS(G(i,j-1)) = RHS(G(i,j-1))-M(topN).*BC.top.c./(BC.top.b/2 + BC.top.a/dy); 

% Bottom boundary
j=1;
i=2:Nx+1;
bottom = reshape(sub2ind(ms, G(i,j+1), G(i,j+1)), Nx,1); % bottom boundary cells
bottomS = reshape(sub2ind(ms, G(i,j+1), G(i,j)), Nx, 1); % south cells to bottom boundary cells
M(bottom) = M(bottom)-((BC.bottom.b/2 + BC.bottom.a/dy)./(BC.bottom.b/2 - BC.bottom.a/dy)).*M(bottomS);
RHS(G(i,j+1)) = RHS(G(i,j+1))-M(bottomS).*BC.bottom.c./(BC.bottom.b/2 - BC.bottom.a/dy);

% Right boundary
i=Nx+2;
j=2:Ny+1;
right = reshape(sub2ind(ms, G(i-1,j), G(i-1,j)), Ny,1); % right boundary cells
rightE = reshape(sub2ind(ms, G(i-1,j), G(i,j)), Ny, 1); % east cells to right boundary cells
M(right) = M(right)-((BC.right.b/2 - BC.right.a/dx)./(BC.right.b/2 + BC.right.a/dx))'.*M(rightE);
RHS(G(i-1,j)) = RHS(G(i-1,j))-M(rightE).*(BC.right.c./(BC.right.b/2 + BC.right.a/dx))';

% Left boundary
i = 1;
j=2:Ny+1;
left = reshape(sub2ind(ms, G(i+1,j), G(i+1,j)), Ny,1); % left boundary cells
leftW = reshape(sub2ind(ms, G(i+1,j), G(i,j)), Ny, 1); % west cells to left boundary cells
M(left) = M(left)-((BC.left.b/2 + BC.left.a/dx)./(BC.left.b/2 - BC.left.a/dx))'.*M(leftW);
RHS(G(i+1,j)) = RHS(G(i+1,j))-M(leftW).*(BC.left.c./(BC.left.b/2 - BC.left.a/dx))';

Mout = M(G(2:end-1,2:end-1), G(2:end-1,2:end-1));
RHSout = RHS(reshape(G(2:end-1,2:end-1),Nx*Ny,1));