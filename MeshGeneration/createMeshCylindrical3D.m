function MS = createMeshCylindrical3D(varargin)
% MeshStructure = createMesh3D(Nr, Ntheta, Nz, Radius, theta, height)
% MeshStructure = createMesh3D(facelocationR, facelocationTheta, facelocationZ)
% builds a uniform 3D mesh:
% Nx is the number of cells in x (horizontal) direction
% Ny is the number of cells in y (vertical) direction
% Nz is the number of cells in z (perpendicular) direction
% Lx is the domain length in x direction
% Ly is the domain length in y direction
% Lz is the domain length in z direction
% 
% SYNOPSIS:
%   MeshStructure = buildMesh3D(Nx, Ny, Nz, Lx, Ly, Lz)
% 
% PARAMETERS:
%   Nx: number of cells in the x direction
%   Lx: domain length in x direction
%   Ny: number of cells in the y direction
%   Ly: domain length in y direction
%   Nz: number of cells in the z direction
%   Lz: domain length in z direction
% 
% RETURNS:
%   MeshStructure.
%                 dimensions=3 (3D problem)
%                 numbering: shows the indexes of cellsn from left to right
%                 and top to bottom and back to front
%                 cellsize: x, y, and z elements of the cell size =[Lx/Nx,
%                 Ly/Ny, Lz/Nz]
%                 cellcenters.x: location of each cell in the x direction
%                 cellcenters.y: location of each cell in the y direction
%                 cellcenters.z: location of each cell in the z direction
%                 facecenters.x: location of interface between cells in the
%                 x direction
%                 facecenters.y: location of interface between cells in the
%                 y direction
%                 facecenters.z: location of interface between cells in the
%                 z direction
%                 numberofcells: [Nx, Ny, Nz]
%                                  
% 
% EXAMPLE:
%   Nx = 2;
%   Lx = 1.0;
%   Ny = 3;
%   Ly = 2.0;
%   Nz = 4;
%   Lz = 3.0;
%   m = buildMesh3D(Nx, Ny, Nz, Lx, Ly, Lz);
%   [X, Y, Z] = ndgrid(m.cellcenters.x, m.cellcenters.y, m.cellcenters.z);
%   [Xf, Yf, Zf] = ndgrid(m.facecenters.x, m.facecenters.y, m.facecenters.z);
%   plot3(X(:), Y(:), Z(:), 'or')
%   hold on;
%   plot3(Xf(:), Yf(:), Zf(:), '+b')
%   legend('cell centers', 'cell corners');
%   
% SEE ALSO:
%     buildMesh1D, buildMesh2D, buildMeshCylindrical1D, ...
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

if nargin==6
  % uniform 1D mesh
  Nx=varargin{1};
  Ny=varargin{2};
  Nz=varargin{3};
  Width=varargin{4};
  Depth=varargin{6};
  Tetta=varargin{5};
  if Tetta>2*pi
      warning('Tetta is higher than 2*pi. It is scaled to 2*pi');
      Tetta = 2*pi;
  end
  MS=createMesh3D(Nx,Ny,Nz,Width, Tetta, Depth);
elseif nargin==3
  % nonuniform 1D mesh
  facelocationX=varargin{1};
  facelocationTheta=varargin{2};
  facelocationZ=varargin{3};
  if facelocationTheta(end)>2*pi
      facelocationTheta = facelocationTheta/facelocationTheta(end)*2.0*pi;
      warning('The domain size adjusted to match a maximum of 2*pi.')
  end
  MS=createMesh3D(facelocationX,facelocationTheta,facelocationZ);
end
MS.dimension=3.2;