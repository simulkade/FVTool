function [M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz] = convectionTvdTerm(MeshStructure, u, phi, FL)
% This function uses the TVD scheme to discretize a
% convection term in the form $\grad (u \phi)$ where u is a face vactor
% It also returns the x, y, x parts of the matrix of coefficient.
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

d = MeshStructure.dimension;
switch d
    case 1
        [M, RHS] = convectionTvdTerm1D(MeshStructure, u, phi, FL);
    case 1.5
        M = convectionTvdTermCylindrical1D(MeshStructure, u, phi, FL);
    case 2
        [M, RHS, Mx, My, RHSx, RHSy] = convectionTvdTerm2D(MeshStructure, u, phi, FL);
    case 2.5
        [M, RHS, Mx, My, RHSx, RHSy] = convectionTvdTermCylindrical2D(MeshStructure, u, phi, FL);
    case 3
        [M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz] = convectionTvdTerm3D(MeshStructure, u, phi, FL);
end
