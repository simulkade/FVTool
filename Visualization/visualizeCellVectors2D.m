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

%{
Copyright (c) 2012, 2013, 2014, 2015 Ali Akbar Eftekhari
All rights reserved.
%}
x = phi_cell.domain.cellcenters.x;
y = phi_cell.domain.cellcenters.y;

quiver(x,y,phi_cell.xvalue', phi_cell.yvalue');
axis equal tight
xlabel('Cell centers [x values]');
ylabel('Cell centers [y values]');
%colorbar
