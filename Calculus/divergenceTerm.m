function [RHSdiv, RHSdivx, RHSdivy, RHSdivz] = divergenceTerm(F)
% This function calculates the divergence of a field using its face
% average value and the vector F, which is a face vector
%
% SYNOPSIS:
%   [RHSdiv, RHSdivx, RHSdivy, RHSdivz] = divergenceTerm(F)
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

RHSdivz=[];

d = F.domain.dimension;
switch d
    case 1
        RHSdiv = divergenceTerm1D(F);
    case 1.5
        RHSdiv = divergenceTermCylindrical1D(F);
    case 1.8
        RHSdiv = divergenceTermSpherical1D(F);
    case 2
        [RHSdiv, RHSdivx, RHSdivy] = divergenceTerm2D(F);
    case 2.5
        [RHSdiv, RHSdivx, RHSdivy] = divergenceTermCylindrical2D(F);
    case 2.8
        [RHSdiv, RHSdivx, RHSdivy] = divergenceTermRadial2D(F);
    case 3
        [RHSdiv, RHSdivx, RHSdivy, RHSdivz] = divergenceTerm3D(F);
    case 3.2
        [RHSdiv, RHSdivx, RHSdivy, RHSdivz] = divergenceTermCylindrical3D(F);
end
