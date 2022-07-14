function r = uminus(p)
    s = struct();
    for idx = 1:numel(p.fields)
        field = p.fields{idx};
        s.(field) = -p.S.(field);
    end
    r = CalculableStruct(s); 
