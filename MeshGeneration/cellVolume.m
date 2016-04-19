function cellvol = cellVolume(meshvar)
% cellvol = cellVolume(meshvar)
% returns the volume of each cell as a cell variable
% SYNOPSIS:
%   cellvol = cellVolume(meshvar)
%
% PARAMETERS:
%   MeshStructure: a mesh structure created by buildMesh* functions
%
% RETURNS:
%   cellvol: a (1D, 2D, or 3D) matrix depending on the mesh size
%
% EXAMPLE:
%   m = createMesh2D(3,4, 1.0, 2.0); % creates a mesh
%   cell_vol=cellVolume(m);
%
% SEE ALSO:
%     createFaceVariable, createBC, buildMesh1D,
%     buildMesh2D, buildMesh3D,
%     buildMeshCylindrical1D, buildMeshCylindrical2D,
%     cellBoundary, combineBC

%{
Copyright (c) 2012, 2013, 2014, 2015, 2016 Ali Akbar Eftekhari
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

% check the size of the variable and the mesh dimension
dim = meshvar.dimension;
BC = createBC(meshvar);
switch dim
    case 1
        c=meshvar.cellsize.x(2:end-1);
    case 1.5
        c=2.0*pi()*meshvar.cellsize.x(2:end-1).*meshvar.cellcenters.x;
    case 2
        c=meshvar.cellsize.x(2:end-1)*meshvar.cellsize.y(2:end-1)';
    case 2.5 % cylindrical
        c=2.0*pi()*meshvar.cellcenters.x.*m.cellsize.x(2:end-1)*m.cellsize.y(2:end-1)';
    case 2.8 % radial
        c=meshvar.cellcenters.x.*m.cellsize.x(2:end-1)*m.cellsize.y(2:end-1)';
    case 3
        error('Not available yet');
    case 3.2
        error('Not available yet');
end
cellvol= createCellVariable(meshvar, c, BC);
