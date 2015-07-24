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

