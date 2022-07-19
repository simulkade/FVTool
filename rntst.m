clear all;
addpath(genpath("."));
rehash
rehash toolboxreset
runtests tests/
% runtests("CalculableStructTest","ProcedureName", ["test_properties", "test_fields",...
% "test_add_field", "test_plus_calculable_structs", "test_from_vec", "test_uminus",...
%     "test_minus", "test_times_calculable_struct", "test_rdivide", "test_copy", "test_sum"])
% runtests("CellTableTest", "ProcedureName", ["test_construction", "test_property_get", ...
%     "test_properties_set", "test_add_field", "test_fieldsCompatible", "test_from_struct",...
%     "test_from_array", "test_plus", "test_uminus", "test_minus", "test_times", "test_rdivide",...
%     "test_sum", "test_mass_mol", "test_get_cv", "test_patch_cv"]);

% runtests({'foo/test_first','foo/test_second'})
% % - or -
% runtests("foo","ProcedureName",["test_first" "test_third"])
% runtests("CalculableStructTest","ProcedureName",["test_fields"])
