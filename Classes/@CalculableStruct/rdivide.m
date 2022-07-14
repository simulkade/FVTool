function r = rdivide(p,q)
    %% Works with celltable, structs of labeled scalars or scalars
    s = struct();
    if isa(p, 'CalculableStruct')&&isa(q, 'CalculableStruct')
        CellTable.fieldsCompatible(p,q);
        for idx = 1:numel(p.fields)
            field = p.fields{idx};
            s.(field) = p.S.(field) ./ q.S.(field);
        end
    elseif isa(p, 'CalculableStruct')&&isa(q, 'CellTable')
        t = table();
        for idx = 1:numel(p.fields)
            field = p.fields{idx};
            t.(field) = p.S.(field) ./ q.T.(field);
        end
        r = CellTable(t);
        return
    elseif isa(p, 'CalculableStruct')&&isa(q, 'struct')
        for idx = 1:numel(p.fields)
            field = p.fields{idx};
            s.(field) = p.S.(field) ./ q.(field);
        end
    elseif isa(p, 'struct') && (isa(q, 'CalculableStruct'))
        CellTable.fieldsCompatible(p,q);
        for idx = 1:numel(q.fields)
            field = q.fields{idx};
            s.(field) = p.(field) ./ q.S.(field);
        end
    elseif isa(p, 'CalculableStruct')
        for idx = 1:numel(p.fields)
            field = p.fields{idx};
            s.(field) = p.S.(field) ./ q;
        end
    else 
        for idx = 1:numel(q.fields)
            field = q.fields{idx};
            s.(field) = p ./ q.S.(field);
        end
    end
    r = CalculableStruct(s);
end
