function visualizeMesh(MS)
%VISUALIZECELLS plots the values of cell variable phi.value
%
% SYNOPSIS:
%   visualizeCells(MS)
%
% PARAMETERS:
%  MS: Mesh Structure
%
% RETURNS:
%    Nothing
%
% EXAMPLE:
%
% SEE ALSO:
%

% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file

d = MS.dimension;
switch d
    case 1
        plot(MS.facecenters.x, zeros(size(MS.facecenters.x)), '-+');
        title('1D Cartesian Grid')
    case 1.5
        L = MS.facecenters.x(end);
        [TH,R] = meshgrid(linspace(-pi/16, pi/16, 3), MS.facecenters.x);
        [X,Y] = pol2cart(TH,R);
        h = polar([0 2*pi], [0 L]);
        delete(h);
        hold on
        pcolor(X,Y,zeros(size(X)));
        title('1D Radial Grid')
        hold off
    case 2
        phi.value(:,1) = 0.5*(phi.value(:,1)+phi.value(:,2));
        phi.value(1,:) = 0.5*(phi.value(1,:)+phi.value(2,:));
        phi.value(:,end) = 0.5*(phi.value(:,end)+phi.value(:,end-1));
        phi.value(end,:) = 0.5*(phi.value(end,:)+phi.value(end-1,:));
        phi.value(1,1) = phi.value(1,2); phi.value(1,end) = phi.value(1,end-1);
        phi.value(end,1) = phi.value(end,2); phi.value(end,end) = phi.value(end,end-1);
        visualizeCells2D(phi);
    case 2.5
        phi.value(:,1) = 0.5*(phi.value(:,1)+phi.value(:,2));
        phi.value(1,:) = 0.5*(phi.value(1,:)+phi.value(2,:));
        phi.value(:,end) = 0.5*(phi.value(:,end)+phi.value(:,end-1));
        phi.value(end,:) = 0.5*(phi.value(end,:)+phi.value(end-1,:));
        phi.value(1,1) = phi.value(1,2); phi.value(1,end) = phi.value(1,end-1);
        phi.value(end,1) = phi.value(end,2); phi.value(end,end) = phi.value(end,end-1);
        visualizeCells2D(phi);
    case 2.8
        [TH,R] = meshgrid(MS.facecenters.y, MS.facecenters.x);
        [X,Y] = pol2cart(TH,R);
        h = polar([0 2*pi], [0 L]);
        delete(h);
        hold on
        pcolor(X,Y,phi.value)
        colorbar
        hold off
    case 3
        phi.value = phi.value(2:end-1,2:end-1,2:end-1);
        visualizeCells3D(phi);
    case 3.2
        phi.value = phi.value(2:end-1,2:end-1,2:end-1);
        visualizeCellsCylindrical3D(phi);
end

end
