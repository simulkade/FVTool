function [BCMatrix, BCRHS] = boundaryCondition(MeshStructure, BC)
% creates the matrix of coefficient and RHS vector
% 
% SYNOPSIS:
%   [BCMatrix, BCRHS] = boundaryCondition(MeshStructure, BC)
% 
% PARAMETERS:
%   MeshStructure  - a mesh structure created by buildMesh* functions
%   BC             - Right hand side values of the boundary condition
%   equations
% 
% RETURNS:
%   BCMatrix  - an square sparse matrix
%   BCRHS     - a column vector of values
% 
% EXAMPLE:
%   L = 1.0; % length of a 1D domain
%   Nx = 5; % Number of grids in x direction
%   m = buildMesh1D(Nx, L);
%   BC = createBC(m); % all Neumann boundaries
%   [Mbc, RHSbc] = boundaryCondition(m, BC);
%   spy(Mbc); % see inside the boundary matrix of coeffiecients
% SEE ALSO:
%     createBC, buildMesh1D, buildMesh2D, buildMesh3D, 
%     buildMeshCylindrical1D, buildMeshCylindrical2D,
%     cellBoundary, combineBC

%{
Copyright (c) 2012, 2013, 2014 Ali Akbar Eftekhari
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

d = MeshStructure.dimension;
if (d ==1) || (d==1.5)
	[BCMatrix, BCRHS] = boundaryCondition1D(MeshStructure, BC);
elseif (d == 2) || (d == 2.5)
	[BCMatrix, BCRHS] = boundaryCondition2D(MeshStructure, BC);
elseif (d == 2.8)
    [BCMatrix, BCRHS] = boundaryConditionRadial2D(MeshStructure, BC);
elseif (d == 3)
    [BCMatrix, BCRHS] = boundaryCondition3D(MeshStructure, BC);
elseif (d == 3.2)
    [BCMatrix, BCRHS] = boundaryConditionCylindrical3D(MeshStructure, BC);
end
