function [M, Mx, My, Mz] = diffusionTerm(D)
% This function uses the central difference scheme to discretize a 2D
% diffusion term in the form \grad . (D \grad \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%   [M, Mx, My, Mz] = diffusionTerm(D)
%
% PARAMETERS:
%   D   - diffusion coefficient, FaceVariable
%
% RETURNS:
%
%
% EXAMPLE:
%
% SEE ALSO:
%


d = D.domain.dimension;
switch d
    case 1
        M = diffusionTerm1D(D);
    case 1.5
        M = diffusionTermCylindrical1D(D);
    case 1.8
        M = diffusionTermSpherical1D(D);
    case 2
        [M, Mx, My] = diffusionTerm2D(D);
    case 2.5
        [M, Mx, My] = diffusionTermCylindrical2D(D);
    case 2.8
        [M, Mx, My] = diffusionTermRadial2D(D);
    case 3
        [M, Mx, My, Mz] = diffusionTerm3D(D);
    case 3.2
        [M, Mx, My, Mz] = diffusionTermCylindrical3D(D);
end
