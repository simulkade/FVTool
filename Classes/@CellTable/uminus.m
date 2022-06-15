function r = uminus(p)
    t = table();
    for idx = 1:numel(p.fields)
        field = p.fields{idx};
        t.(field) = -p.T.(field);
    end
    r = CellTable(t); 
end
    
