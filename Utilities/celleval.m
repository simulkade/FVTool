function r = celleval(f, varargin)
%FUNCEVAL calculates the value of the function f at the cell variable phi
% f must be able to accept matrices
% f: function handle
% phi: input to the function f
% r: output of the function f
% f can handle 6 input phi variables, but they should all be of the same
% class
%
% SYNOPSIS:
%   r = celleval(f, varargin)
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

% % Written by Ali A. Eftekhari
% See the license file

% analyze phi to find the dimension of the matrix

switch nargin
    case 2
        r= CellVariable(varargin{1}.domain, ...
        f(varargin{1}.value));
    case 3
        r= CellVariable(varargin{1}.domain, ...
        f(varargin{1}.value, varargin{2}.value));
    case 4
        r= CellVariable(varargin{1}.domain, ...
        f(varargin{1}.value, varargin{2}.value, varargin{3}.value));
    case 5
        r= CellVariable(varargin{1}.domain, ...
        f(varargin{1}.value, varargin{2}.value, varargin{3}.value, varargin{4}.value));
    case 6
        r= CellVariable(varargin{1}.domain, ...
        f(varargin{1}.value, varargin{2}.value, varargin{3}.value, varargin{4}.value, varargin{5}.value));
end
