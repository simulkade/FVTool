function ct = power(p,q)
    ct = CellTable.from_array(p.mesh, p.A.^q, p.field_struct);
end
