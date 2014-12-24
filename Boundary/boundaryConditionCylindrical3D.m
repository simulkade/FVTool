function [BCMatrix, BCRHS] = boundaryConditionCylindrical3D(MeshStructure, BC)
% It creates the matrix of coefficient based on the BC structureprovided 
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

% Note: I use a for loop here fr more readability of the code!

% extract data from the mesh structure
G = MeshStructure.numbering;
Nxyz = MeshStructure.numberofcells;
Nx = Nxyz(1); Ntetta = Nxyz(2); Nz = Nxyz(3);
d = MeshStructure.cellsize;
dx = d(1); dtetta = d(2); dz = d(3);
rp = repmat(MeshStructure.cellcenters.x', 1, Nz);

% number of boundary nodes (axact number is 2[(m+1)(n+1)*(n+1)*(p+1)+(m+1)*p+1]:
nb = 2*((Nx+1)*(Ntetta+1)+(Nx+1)*(Nz+1)+(Ntetta+1)*(Nz+1));

% define the vectors to be used for the creation of the sparse matrix
ii = zeros(nb,1);
jj = zeros(nb,1);
s = zeros(nb,1);

% define the RHS column vector
BCRHS = zeros((Nx+2)*(Ntetta+2)*(Nz+2), 1);

% assign value to the corner nodes (useless cells)
q = 1:8;
ii(q) = MeshStructure.corner; jj(q) = MeshStructure.corner;
s(q) = 1; BCRHS(MeshStructure.corner) = 0;

% assign values to the edges (useless cells)
q = q(end)+(1:length(MeshStructure.edge));
ii(q) = MeshStructure.edge; jj(q) = MeshStructure.edge;
s(q) = 1; BCRHS(MeshStructure.edge) = 0;

% Assign values to the boundary condition matrix and the RHS vector based
% on the BC structure
if (BC.top.periodic ==0) && (BC.bottom.periodic == 0)
    % top boundary
    j=Ntetta+2;
    i=2:Nx+1;
    k=2:Nz+1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = BC.top.b/2 + BC.top.a./(dtetta*rp);
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j-1,k); s(q) = BC.top.b/2 - BC.top.a./(dtetta*rp);
    BCRHS(G(i,j,k)) = BC.top.c; 

    % Bottom boundary
    j=1;
    i=2:Nx+1;
    k=2:Nz+1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j+1,k);  s(q) = -(BC.bottom.b/2 + BC.bottom.a./(dtetta*rp)); % consider the reverse direction of normal
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k); s(q) = -(BC.bottom.b/2 - BC.bottom.a./(dtetta*rp)); % consider the reverse direction of normal
    BCRHS(G(i,j,k)) = -(BC.bottom.c);
elseif (BC.top.periodic ==1) || (BC.bottom.periodic == 1) % periodic
    % top boundary
    j=Ntetta+2;
    i=2:Nx+1;
    k=2:Nz+1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,2,k); s(q) = -1;
    BCRHS(G(i,j,k)) = 0; 

    % Bottom boundary
    j=1;
    i=2:Nx+1;
    k=2:Nz+1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,Ntetta+1,k); s(q) = -1;
    BCRHS(G(i,j,k)) = 0;
end

if (BC.right.periodic == 0) && (BC.left.periodic == 0)
    % Right boundary
    i=Nx+2;
    j=2:Ntetta+1;
    k=2:Nz+1;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = BC.right.b/2 + BC.right.a/dx;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i-1,j,k); s(q) = BC.right.b/2 - BC.right.a/dx;
    BCRHS(G(i,j,k)) = BC.right.c;

    % Left boundary
    i = 1;
    j=2:Ntetta+1;
    k=2:Nz+1;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i+1,j,k);  s(q) = -(BC.left.b/2 + BC.left.a/dx); % consider the reverse direction of normal
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k); s(q) = -(BC.left.b/2 - BC.left.a/dx); % consider the reverse direction of normal
    BCRHS(G(i,j,k)) = -(BC.left.c);
elseif (BC.right.periodic == 1) || (BC.left.periodic == 1) % periodic
    % Right boundary
    i=Nx+2;
    j=2:Ntetta+1;
    k=2:Nz+1;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(2,j,k); s(q) = -1;
    BCRHS(G(i,j,k)) = 0;

    % Left boundary
    i = 1;
    j=2:Ntetta+1;
    k=2:Nz+1;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(Nx+1,j,k); s(q) = -1;
    BCRHS(G(i,j,k)) = 0;
end

if (BC.front.periodic == 0) && (BC.back.periodic == 0)
    % Back boundary
    k=1;
    i = 2:Nx+1;
    j=2:Ntetta+1;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k+1);  s(q) = -(BC.back.b/2 + BC.back.a/dz); % consider the reverse direction of normal
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k); s(q) = -(BC.back.b/2 - BC.back.a/dz); % consider the reverse direction of normal
    BCRHS(G(i,j,k)) = -(BC.back.c);

    % Front boundary
    k=Nz+2;
    i = 2:Nx+1;
    j=2:Ntetta+1;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = BC.front.b/2 + BC.front.a/dz;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k-1); s(q) = BC.front.b/2 - BC.front.a/dz;
    BCRHS(G(i,j,k)) = BC.front.c;
elseif (BC.front.periodic == 1) || (BC.back.periodic == 1) % periodic
    % Back boundary
    k=1;
    i = 2:Nx+1;
    j=2:Ntetta+1;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,Nz+1); s(q) = -1;
    BCRHS(G(i,j,k)) = 0;

    % Front boundary
    k=Nz+2;
    i = 2:Nx+1;
    j=2:Ntetta+1;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,2); s(q) = -1;
    BCRHS(G(i,j,k)) = 0;
end

% Build the sparse matrix of the boundary conditions
BCMatrix = sparse(ii(1:q(end)), jj(1:q(end)), s(1:q(end)), ...
    (Nx+2)*(Ntetta+2)*(Nz+2), (Nx+2)*(Ntetta+2)*(Nz+2));
