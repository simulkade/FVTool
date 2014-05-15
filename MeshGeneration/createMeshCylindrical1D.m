function MeshStructure = createMeshCylindrical1D(Nr, Lr)
% MeshStructure = buildMeshCylindrical1D(Nr, Lr)
% builds a uniform 1D mesh (1D axial symmetry):
% Nx is the number of cells in r (radial) direction
% Width is the domain radius in r direction
% 
% SYNOPSIS:
%   MeshStructure = buildMeshCylindrical1D(Nr, Lr)
% 
% PARAMETERS:
%   Nr: number of cells in the domain
%   Lr: domain radius
% 
% RETURNS:
%   MeshStructure.
%                 dimensions=1.5 (1D axial symmetry)
%                 numbering: shows the indexes of cellsn from left to right
%                 dx: cell size (=Lr/Nx)
%                 cellcenters.x: location of each cell in the x direction
%                 facecenters.x: location of interface between cells in the
%                 r direction
%                 numberofcells: Nr
%                                  
% 
% EXAMPLE:
%   R = 1.0; % length of the domain
%   Nr = 10; % number of cells in the domain
%   m = buildMeshCylindrical1D(Nr, R);
%   plot(m.cellcenters.x, ones(size(m.cellcenters.x)), 'o', ...
%        m.facecenters.x, ones(size(m.facecenters.x)), '-+');
%   legend('cell centers', 'face centers');
%   
% SEE ALSO:
%     buildMesh1D, buildMesh2D, buildMesh3D, ...
%     buildMeshCylindrical2D, createCellVariable, createFaceVariable

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

% mesh dimension
MeshStructure.dimension = 1.5;

% numbering system of cells, like the single index numbering of Matlab
% +2 is added to account for the ghost cells that are added at 
% the boundaries
MeshStructure.numbering = (1:(Nr+2))'; %numbering system

% cell size is dx
dx = Lr/Nr;

% cell centers position
MeshStructure.cellcenters.x = (1:Nr)*dx-dx/2;

% face centers position
MeshStructure.facecenters.x = (0:Nr)*dx;

% number of cells and cell size in x direction
MeshStructure.numberofcells = Nr;
MeshStructure.cellsize = dx;

