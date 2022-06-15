function s = structProd(svec)
%STRUCTPROD returns the field-by-field products of a structure array
%   S = STRUCTPROD(SARR) takes a structure array and returns a scalar
%   structure with the same fields. The value of each field of S is the
%   product of the values of the corresponding fields of SARR, which must
%   all be numerical scalars.
%
%   See also: prod
fnames = fieldnames(svec);
for k = 1:length(fnames)
    fname = fnames{k};
    s.(fname) = prod([svec.(fname)]);
end
end
