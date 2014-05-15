function [krg, dkrgds] = gaswaterrelperm(s, varargin)
%GASRELPERM returns the relative permeability of the gas phase as a
%function of the gas phase saturation. Other variables can also be
%specified, for instance the critical gas saturation, connate water
%saturaton, etc.
%   s can be gas or liquid saturation and it can be a normal matrix or a
%   face variable structure
% varargin.phase = 'gas' or 'liquid' specifies the phase vor the input s

% check the options value and adapt the saturation variables for the
% calculation
if nargin == 2
    slc = varargin{1}.IrreducibleWaterSaturation;
    sgc = varargin{1}.IrreducibleGasSaturation;
    if strcmpi(varargin{1}.phase, 'gas')
        sg = s;
        dssds = -1/(1-slc-sgc);
    else
        sg = 1 - s;
        dssds = 1/(1-slc-sgc);
    end
else
    sg = s;
    slc = 0;
    sgc = 0;
    dssds = -1/(1-slc-sgc);
end

% calculate the normalized value of the saturation
ss = (1-sg-slc)./(1-slc-sgc);

% load the rel perm functions
kg = @(S)relpermcorey(S);
dkg = @(S)diffrelpermcorey(S, dssds);

% calculate the value of the relperm based on the size of the saturation
% variable
if isstruct(s)==1
    m = size(s.xvalue);
    d = ndims(s.xvalue);
    if (min(m) == 1) % 1D: only x value
        krg.xvalue = kg(ss.xvalue);
        dkrgds.xvalue = dkg(ss.xvalue);
    elseif (d == 2) % 2D: x and y values
        krg.xvalue = kg(ss.xvalue);
        krg.yvalue = kg(ss.yvalue);
        dkrgds.xvalue = dkg(ss.xvalue);
        dkrgds.yvalue = dkg(ss.yvalue);
    else % 3D:x, y, and z values
        krg.xvalue = kg(ss.xvalue);
        krg.yvalue = kg(ss.yvalue);
        krg.zvalue = kg(ss.zvalue);
        dkrgds.xvalue = dkg(ss.xvalue);
        dkrgds.yvalue = dkg(ss.yvalue);
        dkrgds.zvalue = dkg(ss.zvalue);
    end
else
    krg = kg(ss);
    dkrgds = dkg(ss);
end
end

% *************************************************************************
% define the gas relperm function
% change these lines to define a new relperm function for the gas phase,
% or define a new rel perm function and load it in the main function
% Corey type equation; Corey (1954); S is the normalized water saturation
% *************************************************************************
function r = relpermcorey(ss_raw)
% filter the ss values
zone1 = (ss_raw<0); zone2 = (ss_raw>1);
S = ss_raw;
S(zone1) = 0;  S(zone2) = 1;
r = (1-S.^2).*(1-S).^2;
end

function r = diffrelpermcorey(ss_raw,dssds)
% filter the ss values
zone1 = (ss_raw<0); zone2 = (ss_raw>1);
S = ss_raw;
S(zone1) = 0;  S(zone2) = 1;
r = (dssds*(-2*(1-S).*(1-S.^2)-2*(1-S).^2.*S));
r(zone1) =0;    r(zone2) = 0;
end