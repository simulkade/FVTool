function new_cellVar = applyBC(cellVar, BC)
    %% Have to do it as a pure function becaue CellVariable is a value class
    new_cellVar = createCellVariable(cellVar.domain, cellVar.ival, BC);
end

