function options = relpermOptions()
%RELPERMOPIONS creates the options structure for the relperm functions

options.IrreducibleWaterSaturation = 0;
options.IrreducibleGasSaturation = 0;
options.IrreducibleOilSaturation = 0;
options.SortingFactor = 0.457;
options.phase = 'gas';
options.maxOilrelperm = 1;
options.maxWaterrelperm = 1;

end

