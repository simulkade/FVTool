function [M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz] = convectionTvdTerm(u, phi, FL)
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

Mz=[];
d = u.domain.dimension;
switch d
    case 1
        [M, RHS] = convectionTvdTerm1D(u, phi, FL);
    case 1.5
        M = convectionTvdTermCylindrical1D(u, phi, FL);
    case 2
        [M, RHS, Mx, My, RHSx, RHSy] = convectionTvdTerm2D(u, phi, FL);
    case 2.5
        [M, RHS, Mx, My, RHSx, RHSy] = convectionTvdTermCylindrical2D(u, phi, FL);
    case 2.8
        [M, RHS, Mx, My, RHSx, RHSy] = ...
            convectionTvdTermRadial2D(u, phi, FL);
    case 3
        [M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz] = convectionTvdTerm3D(u, phi, FL);
    case 3.2
        [M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz] = ...
            convectionTvdTermCylindrical3D(u, phi, FL);
end
