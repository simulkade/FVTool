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
%}

% domain_length = MeshStructure.numberofcells.*MeshStructure.cellsize;
[X,Y,Z]=meshgrid(phi.domain.cellcenters.y, phi.domain.cellcenters.x, ...
    phi.domain.cellcenters.z);
% Sx = 0.5*domain_length(1);
% Sy = domain_length(2)*[0.25 0.5 0.75];
% Sz = domain_length(3)*[0.25 0.75];
phi.value(1)=phi.value(1)+eps; % to avoid an strange error for assigning color limits
Sx = [phi.domain.cellcenters.x(1) phi.domain.cellcenters.x(end)];
Sy = [phi.domain.cellcenters.y(1) phi.domain.cellcenters.y(end)];
Sz = [phi.domain.cellcenters.z(1) phi.domain.cellcenters.z(end)];
slice(X,Y,Z, phi.value, Sy, Sx, Sz);
xlabel('Cell centers [y vlaues]'); % this is correct [matrix not rotated]
ylabel('Cell centers [x vlaues]'); % this is correct [matrix not rotated]
zlabel('Cell centers [z vlaues]');
axis equal tight
colorbar
