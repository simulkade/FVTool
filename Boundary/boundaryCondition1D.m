function [BCMatrix, BCRHS] = boundaryCondition1D(MeshStructure, BC)
% It creates the matrix of coefficient based on the BC structure provided 
% by the user. It also generates the right hand side vector of the linear
% system of equations
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

% number of boundary nodes:
nb = 2;

% define the vectors to be used for the creation of the sparse matrix
ii = zeros(nb,1);
jj = zeros(nb,1);
s = zeros(nb,1);

% define the RHS column vector
BCRHS = zeros(Nx+2, 1);
q = 0;
% Assign values to the boundary condition matrix and the RHS vector based
% on the BC structure
if (BC.right.periodic == 0) || (BC.left.periodic == 0) % non-periodic boundary condition
    % Right boundary
    i = Nx+2;
    q=q+1; 
    ii(q) = G(i);  jj(q) = G(i);  s(q) = BC.right.b/2 + BC.right.a/dx;
    q=q+1;
    ii(q) = G(i);  jj(q) = G(i-1); s(q) = BC.right.b/2 - BC.right.a/dx;
    BCRHS(G(i)) = BC.right.c;

    % Left boundary
    i = 1;
    q=q+1;
    ii(q) = G(i);  jj(q) = G(i+1);  s(q) = -(BC.left.b/2 + BC.left.a/dx);
    q=q+1;
    ii(q) = G(i);  jj(q) = G(i); s(q) = -(BC.left.b/2 - BC.left.a/dx);
    BCRHS(G(i)) = -(BC.left.c);
elseif (BC.right.periodic == 1) || (BC.left.periodic == 1) % periodic boundary condition
    % Right boundary
    i = Nx+2;
    q=q+1; 
    ii(q) = G(i);  jj(q) = G(i);  s(q) = 1;
    q=q+1;
    ii(q) = G(i);  jj(q) = G(2); s(q) = -1; 
    BCRHS(G(i)) = 0;
    % Left boundary
    i = 1;
    q=q+1;
    ii(q) = G(i);  jj(q) = G(i);  s(q) = 1;
    q=q+1;
    ii(q) = G(i);  jj(q) = G(Nx+1); s(q) = -1;
    BCRHS(G(i)) = 0;
end

% Build the sparse matrix of the boundary conditions
BCMatrix = sparse(ii(1:q), jj(1:q), s(1:q), Nx+2, Nx+2);
