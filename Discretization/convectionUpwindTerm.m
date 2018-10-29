function [M, Mx, My, Mz] = convectionUpwindTerm(u, varargin)
% This function uses the upwind scheme to discretize a 2D
% convection term in the form \grad (u \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%   [M, Mx, My, Mz] = convectionUpwindTerm(u)
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
if nargin>1
    switch d
        case 1
            M = convectionUpwindTerm1D(u, varargin{1});
        case 1.5
            M = convectionUpwindTermCylindrical1D(u, varargin{1});
        case 1.8
            M = convectionUpwindTermSpherical1D(u, varargin{1});
        case 2
            [M, Mx, My] = convectionUpwindTerm2D(u, varargin{1});
        case 2.5
            [M, Mx, My] = convectionUpwindTermCylindrical2D(u, varargin{1});
        case 2.8
            [M, Mx, My] = convectionUpwindTermRadial2D(u, varargin{1});
        case 3
            [M, Mx, My, Mz] = convectionUpwindTerm3D(u, varargin{1});
        case 3.2
            [M, Mx, My, Mz] = convectionUpwindTermCylindrical3D(u, varargin{1});
    end
else
    switch d
        case 1
            M = convectionUpwindTerm1D(u);
        case 1.5
            M = convectionUpwindTermCylindrical1D(u);
        case 1.8
            M = convectionUpwindTermSpherical1D(u);
        case 2
            [M, Mx, My] = convectionUpwindTerm2D(u);
        case 2.5
            [M, Mx, My] = convectionUpwindTermCylindrical2D(u);
        case 2.8
            [M, Mx, My] = convectionUpwindTermRadial2D(u);
        case 3
            [M, Mx, My, Mz] = convectionUpwindTerm3D(u);
        case 3.2
            [M, Mx, My, Mz] = convectionUpwindTermCylindrical3D(u);
    end
end