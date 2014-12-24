function MeshStructure = buildMeshRadial2D(Nr, Ntetta, Lr, Tetta)
% MeshStructure = buildMeshRadial2D(Nr, Ntetta, Lr, Tetta)
% builds a uniform 2D mesh on a Radia coordinate:
% imagine a top view slice of pie 
% Nr is the number of cells in r (radial) direction
% Ntetta is the number of cells in tetta direction
% Lr is the domain length in r direction
% Tetta is the domain length in tetta direction
% 
% SYNOPSIS:
%   MeshStructure = buildMeshRadial2D(Nr, Ntetta, Lr, Tetta)
% 
% PARAMETERS:
%   Nr: number of cells in the x direction
%   Lr: domain length in x direction
%   Ntetta: number of cells in the y direction
%   Tetta: domain length in y direction 0<Tetta<=2*pi
% 
% RETURNS:
%   MeshStructure.
%                 dimensions=2.8 (2D problem, radial coordinate)
%                 numbering: shows the indexes of cellsn from left to right
%                 and top to bottom
%                 cellsize: r and y elements of the cell size =[Lr/Nr,
%                 Ly/Ntetta]
%                 cellcenters.x: location of each cell in the r direction
%                 cellcenters.y: location of each cell in the y direction
%                 facecenters.x: location of interface between cells in the
%                 r direction
%                 facecenters.y: location of interface between cells in the
%                 y direction
%                 numberofcells: [Nr, Ntetta]
%                                  
% 
% EXAMPLE:
%   Nr = 5;
%   Ntetta = 7;
%   R = 10;
%   Ly = 2*pi;
%   m = buildMeshRadial2D(Nx, Ntetta, Lx, Ly);
%   [X, Y] = ndgrid(m.cellcenters.x, m.cellcenters.y);
%   [Xf,Yf]=ndgrid(m.facecenters.x, m.facecenters.y);
%   plot(X, Y, 'or', ...
%        Xf, Yf, '-b', Xf', Yf', '-b');
%   
% SEE ALSO:
%     buildMesh1D, buildMesh3D, buildMeshCylindrical1D, ...
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
AND ANtetta EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
CONTRIBUTORS BE LIABLE FOR ANtetta DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANtetta THEORY OF 
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
NEGLIGENCE OR OTHERWISE) ARISING IN ANtetta WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%}

% mesh dimension
MeshStructure.dimension = 2.8;

% numbering system of cells, like the single index numbering of Matlab
% +2 is added to account for the ghost cells that are added at 
% the boundaries
G = reshape(1:(Nr+2)*(Ntetta+2),Nr+2,Ntetta+2); %numbering system
MeshStructure.numbering = G;

% cell size is dx*dy
dx = Lr/Nr;
if Tetta>2*pi
    warning('Tetta is higher than 2*pi. It is scaled to 2*pi');
    Tetta = 2*pi;
end
dy = Tetta/Ntetta;

% cell centers position
MeshStructure.cellcenters.x = (1:Nr)*dx-dx/2;
MeshStructure.cellcenters.y = (1:Ntetta)*dy-dy/2;

% face centers position
MeshStructure.facecenters.x = (0:Nr)*dx;
MeshStructure.facecenters.y = (0:Ntetta)*dy;

% number of cells and cell size in x and y direction
MeshStructure.numberofcells = [Nr, Ntetta];
MeshStructure.cellsize = [dx, dy];

% boundary indexes
% corner points and edges index
c = G([1 end], [1 end]);
MeshStructure.corner = c(:);