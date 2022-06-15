classdef CellTable < dynamicprops
    % WIP: the idea is to get a data strcture to hold several labeled
    % CellVariables that one may want to treat in the same way and make some
    % algebraic calculations on. Primary use case is concentration of several
    % species.
    % Implementation is an inner table to hold the cellVariables, because matlab
    % does not allow subclassing of table...

    properties
        T  % inner table
    end

    properties (Dependent)
        fields
        mesh
    end

    methods
        function ct = CellTable(inner_table)
            arguments
                inner_table table = table()
            end
            ct.T = inner_table;
            fields = CellTable.safe_fields(ct.T);

            for idx = 1:numel(CellTable.safe_fields(ct.T))
                field = fields{idx};
                ct.add_prop_from_table(field);
            end
        end

        function add_prop_from_table(self, name)
                prop = addprop(self, name);
                prop.Dependent = true;
                prop.GetMethod = @(obj)CellTable.getDynamicProp(obj,name);
                prop.SetMethod = @(obj,val)CellTable.setDynamicProp(obj,name,val);
        end

        function add_field(self, name, cell_variable)
            self.T.(name) = cell_variable;
            self.add_prop_from_table(name);
        end

        function f = get.fields(self)
            f = self.T.Properties.VariableNames;
        end

        function m = get.mesh(self)
            if numel(self.T) == 0
                 m = "Empty CellTable";
            else 
                ex_field = self.fields{1};
                m = self.T.(ex_field).domain;
            end
        end

        function str = repr(self)
            if numel(self.T) == 0
                 str = "Empty CellTable";
            else 
                m = self.mesh;
                dimension = m.dimension;
                if dimension > 1
                    disp('Rep for dim > 1 not implemented');
                    return
                end
                n_cells = m.dims(1);
                x = [-inf; m.cellcenters.x; +inf];
                cell_n = [0; (1:n_cells)'; n_cells+1]; 
                tab = table(cell_n, x);
                for idx = 1:numel(self.fields)
                    field = self.fields{idx};
                    tab.(field) = self.T.(field).value;
                end
                disp(tab);
            end
        end

        function res = sum(self)
            % Does horizontal sum
            res = createCellVariable(self.mesh, 0);
            for idx = 1:numel(self.fields)
                field = self.fields{idx};
                res.value = res.value + self.T.(field).value;
            end
        end

        function self = apply_BC(self, bc)
            sides = {'left', 'right'};
            % Potentially for future: add top, bottom, front, back...
            % when I want to add other domensions and space representations
            for idx = 1:numel(self.fields)
                field = self.fields{idx};
                cell_var = self.T.(field);
                loc_bc = createBC(self.mesh);
                for idx_side = 1:numel(sides)
                    side = sides{idx_side};
                    rob = bc.(side);
                    try
                        a = rob.a.(field);
                    catch
                        a = rob.a;
                    end

                    try
                        b = rob.b.(field);
                    catch
                        b = rob.b;
                    end

                    try
                        c = rob.c.(field);
                    catch
                        c = rob.c;
                    end

                    loc_bc.(side) = robin(a,b,c);
                end

                self.(field) = cell_var.apply_BC(loc_bc);
            end
        end
            
    end

    methods (Static)

        function f = safe_fields(p)
            if isa(p, 'CellTable')
                f = p.fields;
            elseif isa(p, 'table')
                f = p.Properties.VariableNames;
            elseif isa(p, 'struct')
                f = fieldnames(p);
            else
                error("CellTable.safe_fields not implemented for this type");
            end
        end

        function ctablesCompatible(p,q)
            if (numel(setdiff(CellTable.safe_fields(p), ...
                    CellTable.safe_fields(q))) > 0)
                error("Not same set of labels for CellTable operation");
            end
        end

        function val = getDynamicProp(obj,name)
          val = obj.T.(name);
        end

        function setDynamicProp(obj,name,val)
          obj.T.(name) = val;
        end

        function ct = from_struct(mesh, my_struct)
            ct = CellTable();
            var_names = CellTable.safe_fields(my_struct);
                for idx = 1:numel(var_names) 
                    var_name = var_names{idx};
                    ct.add_field(var_name, createCellVariable(mesh, ...
                    my_struct.(var_name)));
                end
        end


    end
end

