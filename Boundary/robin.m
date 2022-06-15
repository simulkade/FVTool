function s = robin(a, b, c)
%% Utility function to create boundary conditions
    s.a = a;
    s.b = b;
    s.c = c;
    s.periodic = 0;
end
