function visualizeCellVectorsRadial2D(phi_cell)
%VISUALIZECELLS plots the values of cell variable phi
%
% SYNOPSIS:
%  visualizeGradsRadial2D(phi)
%
% PARAMETERS:
%  phi_cell: A CellVector variable
%
% RETURNS:
%  None
%
% EXAMPLE:
%
% SEE ALSO:
%

%{
Copyright (c) 2012, 2013, 2014, 2015, Ali Akbar Eftekhari
All rights reserved. Please see the license file.
%}

% phi_cell=gradientCellTerm(phi);
L = phi_cell.domain.cellcenters.x(end);
x = phi_cell.domain.cellcenters.x;
y = phi_cell.domain.cellcenters.y;
[TH,R] = meshgrid(y, x);
[X,Y] = pol2cart(TH,R);
h = polar([0 2*pi], [0 L]);
delete(h);
quiver(X,Y,phi_cell.xvalue.*cos(TH)-phi_cell.yvalue.*sin(TH), phi_cell.xvalue.*sin(TH)+phi_cell.yvalue.*cos(TH))
