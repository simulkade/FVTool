function [pc, dpcds,d2pcds2] = wateroilcapillarypressure(s, varargin)
%CAPILLARYPRESSURE gives the capillary pressure and its derivative with
%respect to the input s: saturation
% default phase for s is the water phase


% check the options value and adapt the saturation variables for the
% calculation
if nargin == 2
    slc = varargin{1}.IrreducibleWaterSaturation;
    soc = varargin{1}.IrreducibleOilSaturation;
    np = varargin{1}.OilWaterExponent;
    pc0 = varargin{1}.OilWaterPc0;
    if strcmpi(varargin{1}.phase, 'water')
        sw = s;
        dssds = -1/(1-slc-soc);
    else
        sw = 1-s;
        dssds = 1/(1-slc-soc);
    end
else
    sw = s;
    slc = 0;
    soc = 0;
    np = 1;
    pc0 = 500; % [Pa]
    dssds = -1/(1-slc-soc);
end

% calculate the normalized value of the saturation
ss = (1-sw-slc)./(1-slc-soc);

% % load the pc functions
p = @(S)cappress(S, pc0, np);
dp = @(S)diffcappress(S, pc0, np, dssds);
d2p = @(S)diff2cappress(S, pc0, np, dssds);

% calculate the value of the relperm based on the size of the saturation
% variable
if isstruct(s)==1
    m = size(s.xvalue);
    d = ndims(s.xvalue);
    if (min(m) == 1) % 1D: only x value
        pc.xvalue = p(ss.xvalue);
        dpcds.xvalue = dp(ss.xvalue);
        d2pcds2.xvalue = d2p(ss.xvalue);
    elseif (d == 2) % 2D: x and y values
        pc.xvalue = p(ss.xvalue);
        pc.yvalue = p(ss.yvalue);
        dpcds.xvalue = dp(ss.xvalue);
        dpcds.yvalue = dp(ss.yvalue);
        d2pcds2.xvalue = d2p(ss.xvalue);
        d2pcds2.yvalue = d2p(ss.yvalue);
    else % 3D:x, y, and z values
        pc.xvalue = p(ss.xvalue);
        pc.yvalue = p(ss.yvalue);
        pc.zvalue = p(ss.zvalue);
        dpcds.xvalue = dp(ss.xvalue);
        dpcds.yvalue = dp(ss.yvalue);
        dpcds.zvalue = dp(ss.zvalue);
        d2pcds2.xvalue = d2p(ss.xvalue);
        d2pcds2.yvalue = d2p(ss.yvalue);
        d2pcds2.zvalue = d2p(ss.zvalue);
    end
else
    pc = p(ss);
    dpcds = dp(ss);
    d2pcds2 = d2p(ss);
end
end


% *************************************************************************
% define the capillary pressure function
% change these two lines to define a new pc function
% *************************************************************************

function r = cappress(ss_raw,  pc0, np)
% filter the ss values
zone1 = (ss_raw<0); zone2 = (ss_raw>1);
S = ss_raw;
S(zone1) = 0;  S(zone2) = 1;
r = pc0*S.^np;
end

function r = diffcappress(ss_raw,  pc0, np, dssds)
% filter the ss values
zone1 = (ss_raw<0); zone2 = (ss_raw>1);
S = ss_raw;
S(zone1) = 0;  S(zone2) = 1;
r = dssds*pc0*np*(S+eps).^(np-1);
% r(zone1) =0;    r(zone2) = 0;
end

function r = diff2cappress(ss_raw,  pc0, np, dssds)
% filter the ss values
zone1 = (ss_raw<0); zone2 = (ss_raw>1);
S = ss_raw;
S(zone1) = 0;  S(zone2) = 1;
r = (dssds^2)*pc0*np*(np-1)*(S+eps).^(np-2);
% r(zone1) =0;    r(zone2) = 0;
end
