function visualizeCellVectorsCylindrical3D(phi_cell)
%VISUALIZECELLS plots the values of cell variable phi
%
% SYNOPSIS:
%   visualizeCellsCylindrical3D(phi)
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

% domain_length = MeshStructure.numberofcells.*MeshStructure.cellsize;
d = phi_cell.domain.dims;
[TH,R,Z]=meshgrid(phi_cell.domain.cellcenters.y, ...
    phi_cell.domain.cellcenters.x, phi_cell.domain.cellcenters.z);
[X,Y,Z]=pol2cart(TH(:),R(:),Z(:)); % in Matlab: [X,Y,Z]=pol2cart(TH,R,Z);
X=reshape(X,d(1),d(2),d(3)); % these three lines work in Matlab and Octave
Y=reshape(Y,d(1),d(2),d(3));
Z=reshape(Z,d(1),d(2),d(3));
quiver3(X,Y,Z,phi_cell.xvalue.*cos(TH)-phi_cell.yvalue.*sin(TH), ...
phi_cell.xvalue.*sin(TH)+phi_cell.yvalue.*cos(TH), phi_cell.zvalue);
axis equal tight
view([60 25]);
% alpha(0.8); % deactivated because it is not implemented in octave
%colorbar;
% xlabel('Cell centers [x vlaues]');
% ylabel('Cell centers [y vlaues]');
%hold off
