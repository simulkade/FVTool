function facevar = createFaceVariable(MeshStructure, faceval)
% this function creates a face variable based on the geometry and mesh
% size. For instance it can be used to create a gravity field.
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

% check the size of the variable and the mesh dimension
d = MeshStructure.dimension;
mn = MeshStructure.numberofcells;

if (d ==1) || (d==1.5)
	facevar.xvalue = faceval.*ones(mn+1, 1);
elseif (d == 2) || (d == 2.5) || (d == 2.8)
	facevar.xvalue = faceval(1).*ones(mn(1)+1, mn(2));
    facevar.yvalue = faceval(2).*ones(mn(1), mn(2)+1);
elseif (d == 3) || (d==3.2)
    facevar.xvalue = faceval(1).*ones(mn(1)+1, mn(2), mn(3));
    facevar.yvalue = faceval(2).*ones(mn(1), mn(2)+1, mn(3));
    facevar.zvalue = faceval(3).*ones(mn(1), mn(2), mn(3)+1);
end