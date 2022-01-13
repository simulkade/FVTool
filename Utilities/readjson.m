% function res = readjson(fname) 
% reads a JSON file and returns the content as an structure
% it is useful in writing parametrized code with nicely formatted input
% files. See the example folder for more some use cases
function res = readjson(fname) 
    fid = fopen(fname); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    res = jsondecode(str);
end