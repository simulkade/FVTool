function phiBC = cellBoundary1D(MeshStructure, BC, phi)
% This function calculates the value of the boundary cells and add them
% to the variable phi of size (1 .. Nx)
% the output includes the boundary is of the size (1..Nx+2)
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
% Nx = MeshStructure.numberofcells;
dx = MeshStructure.cellsize;

% boundary condition (a d\phi/dx + b \phi = c, a column vector of [d a])
% a (phi(i)-phi(i-1))/dx + b (phi(i)+phi(i-1))/2 = c
% phi(i) (a/dx+b/2) + phi(i-1) (-a/dx+b/2) = c
% Right boundary, i=m+2
% phi(i) (a/dx+b/2) = c- phi(i-1) (-a/dx+b/2)
% Left boundary, i=2
%  phi(i-1) (-a/dx+b/2) = c - phi(i) (a/dx+b/2)
% define the new phi
if (BC.left.periodic==0) && (BC.right.periodic==0)
    phiBC = [(BC.left.c-phi(1)* ...
        (BC.left.a/dx+BC.left.b/2))/(-BC.left.a/dx+BC.left.b/2); phi; ...
        (BC.right.c-phi(end)* ...
        (-BC.right.a/dx+BC.right.b/2))/(BC.right.a/dx+BC.right.b/2)];
else
    phiBC = [phi(end); phi; phi(1)];
end