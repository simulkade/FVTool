function BC = createBC3D(meshvar)
% Creates a boundary condition structure from a mesh structure
% for a 3D structured mesh. The default boundary conditions on 
% all boundaries are Neumann;
% The values of each boundary condition are defined as:
% BC.<position>. :
%               a, b, c, where
% a*grad(phi).n + b*phi = c
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

% Extract number of cells from the mesh structure
Nxyz = meshvar.dims;
Nx = Nxyz(1); Ny = Nxyz(2); Nz = Nxyz(3);

% Define the top, bottom, right, and left boundary conditions 
% (default = Neumann, i.e., a = 1, b = 0, c = 0)
top.a = ones(Nx,Nz);
top.b = zeros(Nx,Nz);
top.c = zeros(Nx,Nz);
top.periodic = 0;

bottom.a = ones(Nx,Nz);
bottom.b = zeros(Nx,Nz);
bottom.c = zeros(Nx,Nz);
bottom.periodic = 0;

right.a = ones(Ny,Nz);
right.b = zeros(Ny,Nz);
right.c = zeros(Ny,Nz);
right.periodic = 0;

left.a = ones(Ny,Nz);
left.b = zeros(Ny,Nz);
left.c = zeros(Ny,Nz);
left.periodic = 0;

front.a = ones(Nx,Ny);
front.b = zeros(Nx,Ny);
front.c = zeros(Nx,Ny);
front.periodic = 0;

back.a = ones(Nx,Ny);
back.b = zeros(Nx,Ny);
back.c = zeros(Nx,Ny);
back.periodic = 0;

BC= BoundaryCondition(meshvar, left, right, bottom, top, back, front);