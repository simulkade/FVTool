function visualizeCells(phi)
%VISUALIZECELLS plots the values of cell variable phi.value
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

%{
Copyright (c) 2012, 2013, Ali Akbar Eftekhari
All rights reserved.
%}

d = phi.domain.dimension;
switch d
    case 1
        phi.value= [0.5*(phi.value(1)+phi.value(2)); phi.value(2:end-1); 0.5*(phi.value(end-1)+phi.value(end))];
        visualizeCells1D(phi);
    case 1.5
        phi.value = [0.5*(phi.value(1)+phi.value(2)); phi.value(2:end-1); 0.5*(phi.value(end-1)+phi.value(end))];
        visualizeCells1D(phi);
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
        phi.value(:,1) = 0.5*(phi.value(:,1)+phi.value(:,2));
        phi.value(1,:) = 0.5*(phi.value(1,:)+phi.value(2,:));
        phi.value(:,end) = 0.5*(phi.value(:,end)+phi.value(:,end-1));
        phi.value(end,:) = 0.5*(phi.value(end,:)+phi.value(end-1,:));
        phi.value(1,1) = phi.value(1,2); phi.value(1,end) = phi.value(1,end-1);
        phi.value(end,1) = phi.value(end,2); phi.value(end,end) = phi.value(end,end-1);
        visualizeCellsRadial2D(phi);
    case 3
        phi.value = phi.value(2:end-1,2:end-1,2:end-1);
        visualizeCells3D(phi);
    case 3.2
        phi.value = phi.value(2:end-1,2:end-1,2:end-1);
        visualizeCellsCylindrical3D(phi);
end

end
