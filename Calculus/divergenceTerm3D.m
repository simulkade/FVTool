function [RHSdiv, RHSdivx, RHSdivy, RHSdivz] = divergenceTerm3D(MeshStructure, faceVariable)
% This function calculates the divergence of a field using its face
% average flux vector facevariable, which is a face vector
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

% extract data from the mesh structure
G = MeshStructure.numbering;
Nxyz = MeshStructure.numberofcells;
Nx = Nxyz(1); Ny = Nxyz(2);  Nz = Nxyz(3);
d = MeshStructure.cellsize;
dx = d(1); dy = d(2); dz = d(3);

% define the vector of cell index
row_index = reshape(G(2:Nx+1,2:Ny+1, 2:Nz+1),Nx*Ny*Nz,1); % main diagonal (only internal cells)

% calculate the flux vector in x and y direction
% note: size(Fx) = [1:m+1, 1:n] and size(Fy) = [1:m, 1:n+1]
Fx = faceVariable.xvalue;
Fy = faceVariable.yvalue;
Fz = faceVariable.zvalue;

% reassign the east, west, north, and south flux vectors for the 
% code readability
Fe = Fx(2:Nx+1,:,:);		Fw = Fx(1:Nx,:,:);
Fn = Fy(:,2:Ny+1,:);     Fs = Fy(:,1:Ny,:);
Ff = Fz(:,:,2:Nz+1);     Fb = Fz(:,:,1:Nz);
% compute the divergence
div_x = (Fe - Fw)/dx;
div_y = (Fn - Fs)/dy;
div_z = (Ff - Fb)/dz;

% define the RHS Vector
RHSdiv = zeros((Nx+2)*(Ny+2)*(Nz+2),1);
RHSdivx = zeros((Nx+2)*(Ny+2)*(Nz+2),1);
RHSdivy = zeros((Nx+2)*(Ny+2)*(Nz+2),1);
RHSdivz = zeros((Nx+2)*(Ny+2)*(Nz+2),1);

% assign the values of the RHS vector
RHSdiv(row_index) = reshape(div_x+div_y+div_z,Nx*Ny*Nz,1);
RHSdivx(row_index) = reshape(div_x,Nx*Ny*Nz,1);
RHSdivy(row_index) = reshape(div_y,Nx*Ny*Nz,1);
RHSdivz(row_index) = reshape(div_z,Nx*Ny*Nz,1);

