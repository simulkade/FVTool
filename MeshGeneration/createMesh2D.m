function MS = createMesh2D(varargin)
% MeshStructure = buildMesh2D(Nx, Ny, Lx, Ly)
% MeshStructure = buildMesh2D(facelocationX, facelocationY)
% builds a uniform 2D mesh:
% Nx is the number of cells in x (horizontal) direction
% Ny is the number of cells in y (vertical) direction
% Lx is the domain length in x direction
% Ly is the domain length in y direction
% 
% SYNOPSIS:
%   MeshStructure = buildMesh2D(Nx, Ny, Lx, Ly)
% 
% PARAMETERS:
%   Nx: number of cells in the x direction
%   Lx: domain length in x direction
%   Ny: number of cells in the y direction
%   Ly: domain length in y direction
% 
% RETURNS:
%   MeshStructure.
%                 dimensions=2 (2D problem)
%                 numbering: shows the indexes of cellsn from left to right
%                 and top to bottom
%                 cellsize: x and y elements of the cell size =[Lx/Nx,
%                 Ly/Ny]
%                 cellcenters.x: location of each cell in the x direction
%                 cellcenters.y: location of each cell in the y direction
%                 facecenters.x: location of interface between cells in the
%                 x direction
%                 facecenters.y: location of interface between cells in the
%                 y direction
%                 numberofcells: [Nx, Ny]
%                                  
% 
% EXAMPLE:
%   Nx = 5;
%   Ny = 7;
%   Lx = 10;
%   Ly = 20;
%   m = buildMesh2D(Nx, Ny, Lx, Ly);
%   [X, Y] = ndgrid(m.cellcenters.x, m.cellcenters.y);
%   [Xf,Yf]=ndgrid(m.facecenters.x, m.facecenters.y);
%   plot(X, Y, 'or', ...
%        Xf, Yf, '-b', Xf', Yf', '-b');
%   
% SEE ALSO:
%     buildMesh1D, buildMesh3D, buildMeshCylindrical1D, ...
%     buildMeshCylindrical2D, createCellVariable, createFaceVariable

%{
Copyright (c) 2012, 2013, 2014, Ali Akbar Eftekhari
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

if nargin==4
  % uniform 1D mesh
  Nx=varargin{1};
  Ny=varargin{2};
  Width=varargin{3};
  Height=varargin{4};
  % cell size is dx
  dx = Width/Nx;
  dy = Height/Ny;
  G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
  CellSize.x= dx*ones(Nx+2,1);
  CellSize.y= dy*ones(Ny+2,1);
  CellSize.z= [0.0];
  CellLocation.x= [1:Nx]'*dx-dx/2;
  CellLocation.y= [1:Ny]'*dy-dy/2;
  CellLocation.z= [0.0];
  FaceLocation.x= [0:Nx]'*dx;
  FaceLocation.y= [0:Ny]'*dy;
  FaceLocation.z= [0.0];
elseif nargin==2
  % nonuniform 1D mesh
  facelocationX=varargin{1};
  facelocationY=varargin{2};
  facelocationX=facelocationX(:);
  facelocationY=facelocationY(:);
  Nx = length(facelocationX)-1;
  Ny = length(facelocationY)-1;
  G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
  CellSize.x= [facelocationX(2)-facelocationX(1); ...
    facelocationX(2:end)-facelocationX(1:end-1); ...
    facelocationX(end)-facelocationX(end-1)];
  CellSize.y= [facelocationY(2)-facelocationY(1); ...
    facelocationY(2:end)-facelocationY(1:end-1); ...
    facelocationY(end)-facelocationY(end-1)];
  CellSize.z= [0.0];
  CellLocation.x= 0.5*(facelocationX(2:end)+facelocationX(1:end-1));
  CellLocation.y= 0.5*(facelocationY(2:end)+facelocationY(1:end-1));
  CellLocation.z= [0.0];
  FaceLocation.x= facelocationX;
  FaceLocation.y= facelocationY;
  FaceLocation.z= [0.0];
end
c=G([1,end], [1,end]);
MS=MeshStructure(2, ...
  [Nx,Ny], ...
  CellSize, ...
  CellLocation, ...
  FaceLocation, ...
  c(:), ...
  [1]);