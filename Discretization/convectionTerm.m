function [M, Mx, My, Mz] = convectionTerm(u)
% This function uses the central difference scheme to discretize a 2D
% convection term in the form \grad (u \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
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
    	M = convectionTerm1D(u);
    case 1.5
        M = convectionTermCylindrical1D(u);
    case 1.8
        M = convectionTermSpherical1D(u);
    case 2
        [M, Mx, My] = convectionTerm2D(u);
    case 2.5
        [M, Mx, My] = convectionTermCylindrical2D(u);
    case 2.8
        [M, Mx, My] = convectionTermRadial2D(u);
    case 3
        [M, Mx, My, Mz] = convectionTerm3D(u);
    case 3.2
        [M, Mx, My, Mz] = convectionTermCylindrical3D(u);
end
