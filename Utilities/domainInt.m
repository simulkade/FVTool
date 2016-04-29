function S = domainInt(phi)
%INTEGRATE integrate the value of phi over the domain defined by the mesh
%structure using the trapezoidal rule
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

v=internalCells(cellVolume(phi.domain));
c=internalCells(phi);
S=sum(v(:).*c(:));