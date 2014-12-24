function [RHSdiv, RHSdivx, RHSdivy] = ...
            divergenceTermRadial2D(MeshStructure, faceVariable)
% This function calculates the divergence of a field using its face
% average value and the vector u, which is a face vector
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
Nrz = MeshStructure.numberofcells;
Nr = Nrz(1); Ntetta = Nrz(2);
dxdy = MeshStructure.cellsize;
dx = dxdy(1); dtetta = dxdy(2);
rp = repmat(MeshStructure.cellcenters.x', 1, Ntetta);
rf = repmat(MeshStructure.facecenters.x', 1, Ntetta);

% define the vector of cell index
row_index = reshape(G(2:Nr+1,2:Ntetta+1),Nr*Ntetta,1); % main diagonal (only internal cells)

% calculate the flux vector
% note: size(Fx) = [1:m+1, 1:n] and size(Fy) = [1:m, 1:n+1]
Fx = faceVariable.xvalue;
Fy = faceVariable.yvalue;

% reassign the east, west, north, and south flux vectors for the 
% code readability
Fe = Fx(2:Nr+1,:);		Fw = Fx(1:Nr,:);
Fn = Fy(:,2:Ntetta+1);       Fs = Fy(:,1:Ntetta);
re = rf(2:Nr+1,:);         rw = rf(1:Nr,:);

% compute the divergence
div_x = (re.*Fe - rw.*Fw)./(rp*dx);
div_y = (Fn - Fs)./(dtetta*rp);

% define the RHS Vector
RHSdiv = zeros((Nr+2)*(Ntetta+2),1);
RHSdivx = zeros((Nr+2)*(Ntetta+2),1);
RHSdivy = zeros((Nr+2)*(Ntetta+2),1);

% assign the values of the RHS vector
RHSdiv(row_index) = reshape(div_x+div_y,Nr*Ntetta,1);
RHSdivx(row_index) = reshape(div_x,Nr*Ntetta,1);
RHSdivy(row_index) = reshape(div_y,Nr*Ntetta,1);
