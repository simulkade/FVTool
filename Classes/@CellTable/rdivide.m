function r = rdivide(p,q)
    %% Works with celltable, structs of labeled scalars or scalars
    t = table();
    if isa(p, 'CellTable')&&isa(q, 'CellTable')
        CellTable.fieldsCompatible(p,q);
        for idx = 1:numel(p.fields)
            field = p.fields{idx};
            t.(field) = p.T.(field) ./ q.T.(field);
        end
    elseif isa(p, 'CellTable')&&isa(q, 'CellVariable')
        for idx = 1:numel(p.fields)
            field = p.fields{idx};
            t.(field) = p.T.(field) ./ q;
        end
    elseif isa(p, 'CellVariable')&&isa(q, 'CellTable')
        for idx = 1:numel(p.fields)
            field = p.fields{idx};
            t.(field) = q.T.(field) ./ p;
        end
    elseif isa(p, 'CellTable') && (isa(q, 'struct') || isa(q, 'CalculableStruct'))
        CellTable.fieldsCompatible(p,q);
        for idx = 1:numel(p.fields)
            field = p.fields{idx};
            t.(field) = p.T.(field) / q.(field);
        end
    elseif isa(q, 'CellTable') && (isa(p, 'struct') || isa(p, 'CalculableStruct'))
        CellTable.fieldsCompatible(p,q);
        for idx = 1:numel(q.fields)
            field = q.fields{idx};
            t.(field) = q.T.(field) / p.(field);
        end
    elseif isa(p, 'CellTable')
        for idx = 1:numel(p.fields)
            field = p.fields{idx};
            t.(field) = p.T.(field) / q;
        end
    else
        for idx = 1:numel(q.fields)
            field = q.fields{idx};
            t.(field) = p ./ q.T.(field);
        end
    end
    r = CellTable(t);
end
