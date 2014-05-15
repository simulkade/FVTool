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