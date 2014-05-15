function [M_out, RHS_out] = maskCells(meshstruct, M, RHS, cellIndex, cellValue)
% This function masks the specified cells by giving them a constant value.
% It modifies the matrix of coefficient and the RHS vector
% 
% SYNOPSIS:
%   [M_out, RHS_out] = maskCells(meshstruct, M, RHS, cellIndex, cellValue)
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

% extract some data
d = meshstruct.dimension;
domain_size = meshstruct.numberofcells+2; % 2 is added for the ghost cells
M_size = size(M);
if d ==1 || d==1.5
    i = sub2ind([domain_size 1], cellIndex(:,1));
elseif d==2 || d==2.5
    i = sub2ind(domain_size, cellIndex(:,1), cellIndex(:,2));
elseif d==3
    i = sub2ind(domain_size, cellIndex(:,1), cellIndex(:,2), cellIndex(:,2));
end

% define the new masked matrix of coefficients
M_masked = sparse(i,i,1, M_size(1), M_size(2));
RHS_masked = zeros(length(RHS),1);
RHS_masked(i) = cellValue;

% zero the masked rows in the original matrix
M(i, :) = 0;
RHS(i) = 0;

% add the new masked matrix to the modifed original
M_out = M+M_masked;
RHS_out = RHS + RHS_masked;
