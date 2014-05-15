function options = capillaryPressureOptions()
%RELPERMOPIONS creates the options structure for the relperm functions

    options.IrreducibleWaterSaturation = 0;
    options.IrreducibleOilSaturation = 0;
    options.SortingFactor = 2;
    options.coef = 0.5;
    options.IFT = 0.03; % [N/m]
    options.porosity = 0.3; % [-]
    options.permeability = 1e-12; % [m^2]
    options.phase = 'water';
    options.OilWaterExponent = 1;
    options.OilWaterPc0 = 100; % [Pa]

end

