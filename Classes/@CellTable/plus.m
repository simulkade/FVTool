function r = plus(p,q)
    %% Works with celltable, structs of labeled scalars or scalars
    if isa(p, 'CellTable')&&isa(q, 'CellTable')
        r = CellTable.from_array(p.mesh, p.A + q.A, p.field_struct);
    elseif isa(p, 'CellTable')&&isa(q, 'CellVariable')
        r = CellTable.from_array(p.mesh, p.A + reshape(q.value, [1 1 p.mesh.dims(1)+2]), ...
        p.field_struct);
    elseif isa(p, 'CellVariable')&&isa(q, 'CellTable')
        r = CellTable.from_array(q.mesh, q.A + reshape(p.value, [1 1 q.mesh.dims(1)+2]), ...
        q.field_struct);
    elseif isa(p, 'CellTable') && isa(q, 'CalculableStruct')
        arr_from_v = repmat(reshape(q.V, [1 numel(q.fields)]), [1 1 p.mesh.dims(1)+2]);
        r = CellTable.from_array(p.mesh, p.A + arr_from_v, p.field_struct);
    elseif isa(q, 'CellTable') && isa(p, 'CalculableStruct')
        arr_from_v = repmat(reshape(p.V, [1 numel(p.fields)]), [1 1 q.mesh.dims(1)+2]);
        r = CellTable.from_array(q.mesh, q.A + arr_from_v, q.field_struct);
    elseif isa(p, 'CellTable')
        r = CellTable.from_array(p.mesh, p.A + q, p.field_struct);
    else 
        r = CellTable.from_array(q.mesh, q.A + p, q.field_struct);
    end
end
