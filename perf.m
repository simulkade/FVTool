st1 = CalculableStruct(struct("O2", 1, "N2", 2, "CH4", 1, "H2O", 3, "CO2", 1));
st1.equip_prop();
st2 = CalculableStruct(struct("O2", 1.5, "N2", 0.5, "CH4", 1, "H2O", 2, "CO2", 1));
st2.equip_prop();

N = 10000;

t_add = zeros(1, N);
t_add2 = zeros(1, N);
t_prop = zeros(1, N);

for idx = 1:N   
    tic
    st3 = st1 + st2;
    t_add(idx) = toc;

    tic
    st3.equip_prop();
    t_prop(idx) = toc;

    tic
    st3 = struct();
    for idx2 = 1:numel(st1.fields)
        field = st1.fields(idx2);
        st3.(field) = st1.(field) + st2.(field);
    end
    st3 = CalculableStruct(st3);
    t_add2(idx) = toc;

end
t_add = sum(t_add)
t_add2 = sum(t_add2)
t_prop = sum(t_prop)

    
