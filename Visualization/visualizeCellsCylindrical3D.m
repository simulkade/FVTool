function visualizeCellsCylindrical3D(phi)
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

% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file

% domain_length = MeshStructure.numberofcells.*MeshStructure.cellsize;
d = phi.domain.dims;
[TH,R,Z]=meshgrid(phi.domain.cellcenters.y, ...
    phi.domain.cellcenters.x, phi.domain.cellcenters.z);
[X,Y,Z]=pol2cart(TH(:),R(:),Z(:)); % in Matlab: [X,Y,Z]=pol2cart(TH,R,Z);
X=reshape(X,d(1),d(2),d(3)); % these three lines work in Matlab and Octave
Y=reshape(Y,d(1),d(2),d(3));
Z=reshape(Z,d(1),d(2),d(3));
hold on
% surf three cylinders
surf(squeeze(X(floor(2*d(1)/3),:,:)), squeeze(Y(floor(2*d(1)/3),:,:)), ...
    squeeze(Z(floor(2*d(1)/3),:,:)), squeeze(phi.value(floor(2*d(1)/3),:,:)));
surf(squeeze(X(floor(d(1)/3),:,:)), squeeze(Y(floor(d(1)/3),:,:)), ...
    squeeze(Z(floor(d(1)/3),:,:)), squeeze(phi.value(floor(d(1)/3),:,:)));
surf(squeeze(X(1,:,:)), squeeze(Y(1,:,:)), ...
    squeeze(Z(1,:,:)), squeeze(phi.value(1,:,:)));

% surf three circles
surf(squeeze(X(:,:,1)), squeeze(Y(:,:,1)), ...
    squeeze(Z(:,:,1)), squeeze(phi.value(:,:,1)));
surf(squeeze(X(:,:,floor(d(3)/2))), squeeze(Y(:,:,floor(d(3)/2))), ...
    squeeze(Z(:,:,floor(d(3)/2))), squeeze(phi.value(:,:,floor(d(3)/2))));
surf(squeeze(X(:,:,end)), squeeze(Y(:,:,end)), ...
    squeeze(Z(:,:,end)), squeeze(phi.value(:,:,end)));

% surf two sections
surf(squeeze(X(:,floor(d(2)/3),:)), squeeze(Y(:,floor(d(2)/3),:)), ...
    squeeze(Z(:,floor(d(2)/3),:)), squeeze(phi.value(:,floor(d(2)/3),:)));
surf(squeeze(X(:,floor(2*d(2)/3),:)), squeeze(Y(:,floor(2*d(2)/3),:)), ...
    squeeze(Z(:,floor(2*d(2)/3),:)), squeeze(phi.value(:,floor(2*d(2)/3),:)));

% pcolor(MeshStructure.cellcenters.x, MeshStructure.cellcenters.y, phi')
axis equal tight
view([60 25]);
% alpha(0.8); % deactivated because it is not implemented in octave
colorbar;
% xlabel('Cell centers [x vlaues]');
% ylabel('Cell centers [y vlaues]');
hold off
