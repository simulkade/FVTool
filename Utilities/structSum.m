function s = structSum(mstruct)
    % Returns the sum of all fields of all struct
%%   See also: prod
s = 0;
fnames = fieldnames(mstruct);
for k = 1:length(fnames)
    fname = fnames{k};
    s = s + mstruct.(fname);
end
end
