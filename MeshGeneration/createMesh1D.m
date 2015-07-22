function MS = createMesh1D(varargin)
% MeshStructure = createMesh1D(Nx, Width)
% MeshStructure = createMesh1D(facelocationX)
% builds a uniform 1D mesh:
% Nx is the number of cells in x (horizontal) direction
% Width is the domain length in x direction
% 
% SYNOPSIS:
%   MeshStructure = buildMesh1D(Nx, Width)
% 
% PARAMETERS:
%   Nx: number of cells in the domain
%   Width: domain length
% 
% RETURNS:
%   MeshStructure.
%                 dimensions=1 (1D problem)
%                 numbering: shows the indexes of cellsn from left to right
%                 dx: cell size (=Width/Nx)
%                 cellcenters.x: location of each cell in the x direction
%                 facecenters.x: location of interface between cells in the
%                 x direction
%                 numberofcells: Nx
%                                  
% 
% EXAMPLE:
%   L = 1.0; % length of the domain
%   Nx = 10; % number of cells in the domain
%   m = buildMesh1D(Nx, L);
%   plot(m.cellcenters.x, ones(size(m.cellcenters.x)), 'o', ...
%        m.facecenters.x, ones(size(m.facecenters.x)), '-+');
%   legend('cell centers', 'face centers');
%   
% SEE ALSO:
%     buildMesh2D, buildMesh3D, buildMeshCylindrical1D, ...
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

if nargin==2
  % uniform 1D mesh
  Nx=varargin{1};
  Width=varargin{2};
  % cell size is dx
  dx = Width/Nx;
  CellSize.x= dx*ones(Nx+2,1);
  CellSize.y= [0.0];
  CellSize.z= [0.0];
  CellLocation.x= [1:Nx]'*dx-dx/2;
  CellLocation.y= [0.0];
  CellLocation.z= [0.0];
  FaceLocation.x= [0:Nx]'*dx;
  FaceLocation.y= [0.0];
  FaceLocation.z= [0.0];
elseif nargin==1
  % nonuniform 1D mesh
  facelocationX=varargin{1};
  n=size(facelocationX);
  if n(1)==1
      facelocationX=facelocationX';
  end
  Nx = length(facelocationX)-1;
  CellSize.x= [facelocationX(2)-facelocationX(1); ...
    facelocationX(2:end)-facelocationX(1:end-1); ...
    facelocationX(end)-facelocationX(end-1)];
  CellSize.y= [0.0];
  CellSize.z= [0.0];
  CellLocation.x= 0.5*(facelocationX(2:end)+facelocationX(1:end-1));
  CellLocation.y= [0.0];
  CellLocation.z= [0.0];
  FaceLocation.x= facelocationX;
  FaceLocation.y= [0.0];
  FaceLocation.z= [0.0];
end

MS=MeshStructure(1, ...
  [Nx,1], ...
  CellSize, ...
  CellLocation, ...
  FaceLocation, ...
  [1], ...
  [1]);
