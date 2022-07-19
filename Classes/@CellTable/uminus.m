function r = uminus(p)
    r = CellTable.from_array(p.mesh, -p.A, p.field_struct);
end
    
