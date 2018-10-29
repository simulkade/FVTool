function MS = createMeshRadial2D(varargin)
% MeshStructure = createMeshRadial2D(Nr, Ntetta, Lr, Tetta)
% MeshStructure = createMeshRadial2D(facelocationR, facelocationTetta)
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
%   m = createMeshRadial2D(Nx, Ntetta, Lx, Ly);
%   [X, Y] = ndgrid(m.cellcenters.x, m.cellcenters.y);
%   [Xf,Yf]=ndgrid(m.facecenters.x, m.facecenters.y);
%   plot(X, Y, 'or', ...
%        Xf, Yf, '-b', Xf', Yf', '-b');
%
% SEE ALSO:
%     buildMesh1D, buildMesh3D, buildMeshCylindrical1D, ...
%     buildMeshCylindrical2D, createCellVariable, createFaceVariable

% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file

if nargin==4
  % uniform 1D mesh
  Nx=varargin{1};
  Ny=varargin{2};
  Width=varargin{3};
  Tetta=varargin{4};
  if Tetta>2*pi
    warning('Tetta is higher than 2*pi. It is scaled to 2*pi');
    Tetta = 2*pi;
  end
  MS=createMesh2D(Nx, Ny, Width, Tetta);
elseif nargin==2
  % nonuniform 1D mesh
  facelocationX=varargin{1};
  facelocationY=varargin{2};
  if facelocationY(end)>2*pi
      facelocationY = facelocationY/facelocationY(end)*2.0*pi;
      warning('The domain size adjusted to match a maximum of 2*pi.')
  end
  MS=createMesh2D(facelocationX, facelocationY);
end
MS.dimension=2.8;
