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

% Written by Ali A. Eftekhari
% See the license file

x = [phi.domain.facecenters.x(1); phi.domain.cellcenters.x; phi.domain.facecenters.x(end)];
y = [phi.domain.facecenters.y(1); phi.domain.cellcenters.y; phi.domain.facecenters.y(end)];

pcolor(x, y, phi.value')
axis equal tight
xlabel('Cell centers [x vlaues]');
ylabel('Cell centers [y vlaues]');
colorbar
