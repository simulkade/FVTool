function visualizeCells(MeshStructure, phi)
%VISUALIZECELLS plots the values of cell variable phi
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

d = MeshStructure.dimension;
domain_size = prod(MeshStructure.numberofcells);
phi_size = numel(phi);
switch d
    case 1
        if domain_size==phi_size
            visualizeCells1D(MeshStructure, phi);
        else
            visualizeCells1D(MeshStructure, [0.5*(phi(1)+phi(2)); phi(2:end-1); 0.5*(phi(end-1)+phi(end))]);
        end
    case 1.5
        if domain_size==phi_size
            visualizeCells1D(MeshStructure, phi);
        else            
            visualizeCells1D(MeshStructure, [0.5*(phi(1)+phi(2)); phi(2:end-1); 0.5*(phi(end-1)+phi(end))]);
        end
    case 2
        if domain_size==phi_size
            visualizeCells2D(MeshStructure, phi);
        else
            phi(:,1) = 0.5*(phi(:,1)+phi(:,2));
            phi(1,:) = 0.5*(phi(1,:)+phi(2,:));
            phi(:,end) = 0.5*(phi(:,end)+phi(:,end-1));
            phi(end,:) = 0.5*(phi(end,:)+phi(end-1,:));
            phi(1,1) = phi(1,2); phi(1,end) = phi(1,end-1);
            phi(end,1) = phi(end,2); phi(end,end) = phi(end,end-1);
            visualizeCells2D(MeshStructure, phi);
        end
    case 2.5
        if domain_size==phi_size
            visualizeCells2D(MeshStructure, phi);
        else
            phi(:,1) = 0.5*(phi(:,1)+phi(:,2));
            phi(1,:) = 0.5*(phi(1,:)+phi(2,:));
            phi(:,end) = 0.5*(phi(:,end)+phi(:,end-1));
            phi(end,:) = 0.5*(phi(end,:)+phi(end-1,:));
            phi(1,1) = phi(1,2); phi(1,end) = phi(1,end-1);
            phi(end,1) = phi(end,2); phi(end,end) = phi(end,end-1);
            visualizeCells2D(MeshStructure, phi);
        end
    case 2.8
        if domain_size==phi_size
            visualizeCellsRadial2D(MeshStructure, phi);
        else
            phi(:,1) = 0.5*(phi(:,1)+phi(:,2));
            phi(1,:) = 0.5*(phi(1,:)+phi(2,:));
            phi(:,end) = 0.5*(phi(:,end)+phi(:,end-1));
            phi(end,:) = 0.5*(phi(end,:)+phi(end-1,:));
            phi(1,1) = phi(1,2); phi(1,end) = phi(1,end-1);
            phi(end,1) = phi(end,2); phi(end,end) = phi(end,end-1);
            visualizeCellsRadial2D(MeshStructure, phi);
        end
    case 3
        if domain_size==phi_size
            visualizeCells3D(MeshStructure, phi);
        else
            visualizeCells3D(MeshStructure, phi(2:end-1,2:end-1,2:end-1));
        end
    case 3.2
        if domain_size==phi_size
            visualizeCellsCylindrical3D(MeshStructure, phi);
        else
            visualizeCellsCylindrical3D(MeshStructure, phi(2:end-1,2:end-1,2:end-1));
        end
end

end

