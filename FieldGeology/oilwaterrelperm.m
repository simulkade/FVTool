function [krl, dkrlds] = oilwaterrelperm(s, varargin)
%OILRELPERM returns the relative permeability of the oil phase as a
%function of the oil phase saturation. Other variables can also be
%specified, for instance the critical gas saturation, connate water
%saturaton, etc.
%   s can be gas or liquid saturation and it can be a normal matrix or a
%   face variable structure; default: s is liquid phase saturation
% varargin.phase = 'gas' or 'liquid' specifies the phase for the input s

% check the options value and adapt the saturation variables for the
% calculation
if nargin == 2
    swc = varargin{1}.IrreducibleWaterSaturation;
    soc = varargin{1}.IrreducibleOilSaturation;
    krw0 = varargin{1}.maxWaterrelperm;
    if strcmpi(varargin{1}.phase, 'water')
        so = 1-s;
        dssds = 1/(1-swc-soc);
    else
        so = s;
        dssds = -1/(1-swc-soc);
    end
else
    so = 1-s;
    swc = 0;
    soc = 0;
    krw0 = 1;
    dssds = -1/(1-swc-soc);
end

% calculate the normalized value of the saturation
ss = (1-so-swc)./(1-swc-soc);

% load the rel perm functions
kl = @(S)relpermcorey(S, krw0);
dkl = @(S)diffrelpermcorey(S, krw0, dssds);

% calculate the value of the relperm based on the size of the saturation
% variable
if isstruct(s)==1
    m = size(s.xvalue);
    d = ndims(s.xvalue);
    if (min(m) == 1) % 1D: only x value
        krl.xvalue = kl(ss.xvalue);
        dkrlds.xvalue = dkl(ss.xvalue);
    elseif (d == 2) % 2D: x and y values
        krl.xvalue = kl(ss.xvalue);
        krl.yvalue = kl(ss.yvalue);
        dkrlds.xvalue = dkl(ss.xvalue);
        dkrlds.yvalue = dkl(ss.yvalue);
    else % 3D:x, y, and z values
        krl.xvalue = kl(ss.xvalue);
        krl.yvalue = kl(ss.yvalue);
        krl.zvalue = kl(ss.zvalue);
        dkrlds.xvalue = dkl(ss.xvalue);
        dkrlds.yvalue = dkl(ss.yvalue);
        dkrlds.zvalue = dkl(ss.zvalue);
    end
else
    krl = kl(ss);
    dkrlds = dkl(ss);
end
end

% *************************************************************************
% define the oil relperm function
% change these lines to define a new relperm function for the oil phase
% *************************************************************************

function r = relpermcorey(ss_raw, krw0)
% filter the ss values
zone1 = (ss_raw<0); zone2 = (ss_raw>1);
S = ss_raw;
S(zone1) = 0;  S(zone2) = 1;
r = krw0*(1-S).^2;
end

function r = diffrelpermcorey(ss_raw,krw0,dssds)
% filter the ss values
zone1 = (ss_raw<0); zone2 = (ss_raw>1);
S = ss_raw;
S(zone1) = 0;  S(zone2) = 1;
r = -2*krw0*dssds*(1-S);
% r(zone1) =0;    r(zone2) = 0;
end