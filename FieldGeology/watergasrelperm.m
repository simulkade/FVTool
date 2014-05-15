function [krl, dkrlds] = watergasrelperm(s, varargin)
%LIQUIDRELPERM returns the relative permeability of the water phase as a
%function of the water phase saturation. Other variables can also be
%specified, for instance the critical gas saturation, connate water
%saturaton, etc.
%   s can be gas or liquid saturation and it can be a normal matrix or a
%   face variable structure; default: s is liquid phase saturation
% varargin.phase = 'gas' or 'liquid' specifies the phase for the input s

% check the options value and adapt the saturation variables for the
% calculation
if nargin == 2
    slc = varargin{1}.IrreducibleWaterSaturation;
    sgc = varargin{1}.IrreducibleGasSaturation;
    L = varargin{1}.SortingFactor;
    if strcmpi(varargin{1}.phase, 'gas')
        sg = s;
        dssds = -1/(1-slc-sgc);
    else
        sg = 1-s;
        dssds = 1/(1-slc-sgc);
    end
else
    sg = 1-s;
    slc = 0;
    sgc = 0;
    L = 0.457;
    dssds = 1/(1-slc-sgc);
end

% calculate the normalized value of the saturation
ss = (1-sg-slc)./(1-slc-sgc);

% load the rel perm functions
kl = @(S)relpermvanGenuchten(S, L);
dkl = @(S)diffrelpermvanGenuchten(S, L, dssds);

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
% define the water relperm function
% change these lines to define a new relperm function for the gas phase
% van Genuchten (1980)
% *************************************************************************
function r = relpermvanGenuchten(ss_raw, L)
% filter the ss values
zone1 = (ss_raw<0); zone2 = (ss_raw>1);
s = ss_raw;
s(zone1) = 0;  s(zone2) = 1;
r = (s.^0.5).*(1-(1-s.^(1/L)).^L).^2;
end

function r = diffrelpermvanGenuchten(ss_raw, L, dssds)
% filter the ss values
zone1 = (ss_raw<=0); zone2 = (ss_raw>=1);
s = ss_raw;
s(zone1) = 0;  s(zone2) = 1;
r = dssds*(2*(s+eps).^(1/L-0.5).*(1-s.^(1/L)+eps).^(L-1).*(1-(1-s.^(1/L)).^L)+ ...
    (0.5*(1-(1-s.^(1/L)).^L).^2)./(s+eps).^0.5);
r(zone1) =0;    r(zone2) = 0;
end