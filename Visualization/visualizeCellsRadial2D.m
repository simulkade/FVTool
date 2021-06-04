function visualizeCellsRadial2D(phi)
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

% Written by Ali A. Eftekhari
% See the license file

L = phi.domain.cellcenters.x(end);
x = [phi.domain.facecenters.x(1); phi.domain.cellcenters.x; phi.domain.facecenters.x(end)];
y = [phi.domain.facecenters.y(1); phi.domain.cellcenters.y; phi.domain.facecenters.y(end)];
[TH,R] = meshgrid(y, x);
[X,Y] = pol2cart(TH,R);
h = polar([0 2*pi], [0 L]);
delete(h);
hold on
pcolor(X,Y,phi.value)
colorbar
hold off
