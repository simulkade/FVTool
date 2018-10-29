function r = uminus(p)
%UMINUS: this function applies a unary minus to the x, y, and z values of the structures that I use in
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

r=p;
r.xvalue = -p.xvalue;
r.yvalue = -p.yvalue;
r.zvalue = -p.zvalue;
