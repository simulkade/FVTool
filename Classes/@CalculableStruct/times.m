function r = times(p,q)
    if isa(p, 'CalculableStruct')&&isa(q, 'CalculableStruct')
        r = CalculableStruct.from_vec(p.V .* q.V, p.field_struct);
    elseif isa(p, 'CalculableStruct')&&isa(q, 'CellTable')
        r = q .* p;
    elseif isa(p, 'CalculableStruct')&&isa(q, 'CellVariable')
        arr_from_cs = repmat(reshape(p.V, [1 numel(p.fields)]), [1 1 q.domain.dims(1)+2]);
        arr_from_cv = reshape(q.value, [1 1 q.domain.dims(1)+2]);
        r = CellTable.from_array(q.domain, arr_from_cs .* arr_from_cv, p.field_struct);
    elseif isa(p, 'CellVariable')&&isa(q, 'CalculableStruct')
        r = q .* p;
    elseif isa(p, 'CalculablStruct')
        r = CalculableStruct.from_vec(p.V .* q, p.field_struct);
    else 
        r = CalculableStruct.from_vec(q.V .* p, q.field_struct);
    end
end
