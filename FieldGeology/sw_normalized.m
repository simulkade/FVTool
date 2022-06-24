function sws = sw_normalized(sw, swc, sor)
%SW_NORMALIZED maps the saturation from [0, 1] to [swc, 1-sor]
sws=((sw>swc).*(sw<1-sor).*(sw-swc)./(1-sor-swc)+(sw>=1-sor));
end

