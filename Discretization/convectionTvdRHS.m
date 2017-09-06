function [RHS, RHSx, RHSy, RHSz] = convectionTvdRHS(u, phi, FL)
% This function uses the TVD scheme to discretize a
% convection term in the form $\grad (u \phi)$ where u is a face vactor
% It also returns the x, y, x parts of the matrix of coefficient.
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

% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file

d = u.domain.dimension;
switch d
    case 1
        RHS = convectionTvdRHS1D(u, phi, FL);
    case 1.5
        RHS = convectionTvdRHSCylindrical1D(u, phi, FL);
    case 2
        [RHS, RHSx, RHSy] = convectionTvdRHS2D(u, phi, FL);
    case 2.5
        [RHS, RHSx, RHSy] = convectionTvdRHSCylindrical2D(u, phi, FL);
    case 2.8
        [RHS, RHSx, RHSy] = ...
            convectionTvdRHSRadial2D(u, phi, FL);
    case 3
        [RHS, RHSx, RHSy, RHSz] = convectionTvdRHS3D(u, phi, FL);
    case 3.2
        [RHS, RHSx, RHSy, RHSz] = ...
            convectionTvdRHSCylindrical3D(u, phi, FL);
end
