function r = times(p,q)
    if isa(p, 'CellTable')&&isa(q, 'CellTable')
        r = CellTable.from_array(p.mesh, p.A .* q.A, p.field_struct);
    elseif isa(p, 'CellTable')&&isa(q, 'CellVariable')
        r = CellTable.from_array(p.mesh, p.A .* reshape(q.value, [1 1 p.mesh.dims(1)+2]), ...
        p.field_struct);
    elseif isa(p, 'CellVariable')&&isa(q, 'CellTable')
        r = CellTable.from_array(q.mesh, q.A .* reshape(p.value, [1 1 q.mesh.dims(1)+2]), ...
        q.field_struct);
    elseif isa(p, 'CellTable') && isa(q, 'CalculableStruct')
        arr_from_v = repmat(reshape(q.V, [1 numel(q.fields)]), [1 1 p.mesh.dims(1)+2]);
        r = CellTable.from_array(p.mesh, p.A .* arr_from_v, p.field_struct);
    elseif isa(p, 'CalculableStruct') && isa(q, 'CellTable')
        r = q .* p;
    elseif isa(p, 'CellTable')
        r = CellTable.from_array(p.mesh, p.A .* q, p.field_struct);
    else
        r = CellTable.from_array(q.mesh, q.A .* p, q.field_struct);
    end
end
