function s = robin(a, b, c)
%% Utility function to create boundary conditions
    if class(a) == "struct"
        fields = fieldnames(a); 
        for idx = 1:numel(fields)
            field = fields{idx};
            s.a.(field) = a.(field);
        end
    else
        s.a = a;
    end

    if class(b) == "struct"
        fields = fieldnames(c); 
        for idx = 1:numel(fields)
            field = fields{idx};
            s.b.(field) = b.(field);
        end
    else
        s.b = b;
    end

    if class(c) == "struct"
        fields = fieldnames(c); 
        for idx = 1:numel(fields)
            field = fields{idx};
            s.c.(field) = c.(field);
        end
    else
        s.c = c;
    end

    s.periodic = 0;

end
