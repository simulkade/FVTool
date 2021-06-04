function visualizeCells1D(phi)
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

x = [phi.domain.facecenters.x(1); phi.domain.cellcenters.x; phi.domain.facecenters.x(end)];
plot(x, phi.value)
xlabel('Cell centers [x vlaues]');
ylabel('Cell values');
