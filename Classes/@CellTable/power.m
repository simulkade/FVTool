function ct = power(p,q)

    t = table();
    if isa(p, 'CellTable')&&isa(q, 'double')
        for idx = 1:numel(p.fields)
            field = p.fields{idx};
            t.(field) = p.T.(field).^q;
        end
    else
        error('CellTable: power with type than table not Implementated') 
    end
    ct = CellTable(t);
end
