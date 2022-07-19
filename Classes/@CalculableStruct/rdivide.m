function r = rdivide(p,q)
    %% Works with celltable, structs of labeled scalars or scalars
    if isa(p, 'CalculableStruct')&&isa(q, 'CalculableStruct')
        r = CalculableStruct.from_vec(p.V ./ q.V, p.field_struct);
    elseif isa(p, 'CalculableStruct')&&isa(q, 'CellVariable')
        arr_from_cs = repmat(reshape(p.V, [1 numel(p.fields)]), [1 1 q.domain.dims(1)+2]);
        arr_from_cv = reshape(q.value, [1 1 q.domain.dims(1)+2]);
        r = CellTable.from_array(q.domain, arr_from_cs ./ arr_from_cv, p.field_struct);
    elseif isa(p, 'CellVariable')&&isa(q, 'CalculableStruct')
        arr_from_cs = repmat(reshape(q.V, [1 numel(q.fields)]), [1 1 p.domain.dims(1)+2]);
        arr_from_cv = reshape(p.value, [1 1 p.domain.dims(1)+2]);
        r = CellTable.from_array(p.domain, arr_from_cv ./ arr_from_cs, q.field_struct);
    elseif isa(p, 'CalculableStruct')
        r = CalculableStruct.from_vec(p.V ./ q, p.field_struct);
    else 
        r = CalculableStruct.from_vec(p ./ q.V, q.field_struct);
    end
end
