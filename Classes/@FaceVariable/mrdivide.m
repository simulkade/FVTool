function r = mrdivide(p,q)
%RDIVIDE this function divides the x, y, and z values of the structures that I use in
% the FVtool.
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

if isa(p, 'FaceVariable') && (~isa(q, 'FaceVariable'))
    r = (1/q)*p;
else
    error('FVMtool: Wrong use of mrdivide for a face variable. Try using ./ instead.');
end
