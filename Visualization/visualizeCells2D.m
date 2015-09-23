function visualizeCells2D(phi)
%VISUALIZECELLS plots the values of cell variable phi
%
% SYNOPSIS:
%   visualizeCellVectors2D(phi_cell)
%
% PARAMETERS:
%   phi_cell: CellVector
%
% RETURNS:
%   None
%
% EXAMPLE:
%
% SEE ALSO:
%

%{
Copyright (c) 2012, 2013, 2014, 2015 Ali Akbar Eftekhari
All rights reserved.
%}
x = [phi.domain.facecenters.x(1); phi.domain.cellcenters.x; phi.domain.facecenters.x(end)];
y = [phi.domain.facecenters.y(1); phi.domain.cellcenters.y; phi.domain.facecenters.y(end)];

pcolor(x, y, phi.value')
axis equal tight
xlabel('Cell centers [x vlaues]');
ylabel('Cell centers [y vlaues]');
colorbar
