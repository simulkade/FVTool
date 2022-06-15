function s = structScalarProd(mstruct, scalar)
fnames = fieldnames(mstruct);
for k = 1:length(fnames)
    fname = fnames{k};
    s.(fname) = mstruct.(fname) * scalar;
end
end
