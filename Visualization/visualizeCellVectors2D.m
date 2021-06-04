function visualizeCellVectors2D(phi_cell)
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

x = phi_cell.domain.cellcenters.x;
y = phi_cell.domain.cellcenters.y;

quiver(x,y,phi_cell.xvalue', phi_cell.yvalue');
axis equal tight
xlabel('Cell centers [x values]');
ylabel('Cell centers [y values]');
%colorbar
