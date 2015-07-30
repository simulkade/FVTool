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

%{
Copyright (c) 2012, 2013, Ali Akbar Eftekhari
All rights reserved.

Redistribution and use in source and binary forms, with or 
without modification, are permitted provided that the following 
conditions are met:

    *   Redistributions of source code must retain the above copyright notice, 
        this list of conditions and the following disclaimer.
    *   Redistributions in binary form must reproduce the above 
        copyright notice, this list of conditions and the following 
        disclaimer in the documentation and/or other materials provided 
        with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%}

d = F.domain.dimension;
switch d
    case 1
        RHSdiv = divergenceTerm1D(F);
    case 1.5
        RHSdiv = divergenceTermCylindrical1D(F);
    case 2
        [RHSdiv, RHSdivx, RHSdivy] = divergenceTerm2D(F);
    case 2.5
        [RHSdiv, RHSdivx, RHSdivy] = divergenceTermCylindrical2D(F);
    case 2.8
        [RHSdiv, RHSdivx, RHSdivy] = ...
            divergenceTermRadial2D(F);
    case 3
        [RHSdiv, RHSdivx, RHSdivy, RHSdivz] = divergenceTerm3D(F);
    case 3.2
        [RHSdiv, RHSdivx, RHSdivy, RHSdivz] = divergenceTermCylindrical3D(F);
end
