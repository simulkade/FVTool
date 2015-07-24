function visualizeCells3D(phi)
%VISUALIZECELLS plots the values of cell variable phi
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

% domain_length = MeshStructure.numberofcells.*MeshStructure.cellsize;
[X,Y,Z]=meshgrid(phi.domain.cellcenters.y, phi.domain.cellcenters.x, ...
    phi.domain.cellcenters.z);
% Sx = 0.5*domain_length(1);
% Sy = domain_length(2)*[0.25 0.5 0.75];
% Sz = domain_length(3)*[0.25 0.75];
Sx = [phi.domain.cellcenters.x(1) phi.domain.cellcenters.x(end)];
Sy = [phi.domain.cellcenters.y(1) phi.domain.cellcenters.y(end)];
Sz = [phi.domain.cellcenters.z(1) phi.domain.cellcenters.z(end)];
slice(X,Y,Z, phi.value, Sy, Sx, Sz);
xlabel('Cell centers [y vlaues]'); % this is correct [matrix not rotated]
ylabel('Cell centers [x vlaues]'); % this is correct [matrix not rotated]
zlabel('Cell centers [z vlaues]');
axis equal tight
colorbar
