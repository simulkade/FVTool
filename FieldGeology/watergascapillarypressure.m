function [pc, dpcds] = watergascapillarypressure(s, varargin)
%CAPILLARYPRESSURE gives the capillary pressure and its derivative with
%respect to the input s: saturation
% default phase for s is the water phase


% check the options value and adapt the saturation variables for the
% calculation
if nargin == 2
    slc = varargin{1}.IrreducibleWaterSaturation;
    sgc = varargin{1}.IrreducibleGasSaturation;
    L = varargin{1}.SortingFactor;
    g = varargin{1}.coef;
    sigma = varargin{1}.IFT; % [N/m]
    phi = varargin{1}.porosity; % [-]
    k = varargin{1}.permeability; % [m^2]
    if strcmpi(varargin{1}.phase, 'gas')
        sg = s;
        dssds = -1/(1-slc);
    else
        sg = 1-s;
        dssds = 1/(1-slc);
    end
else
    sg = 1-s;
    slc = 0;
    sgc = 0;
    L = 2;
    g = 0.5;
    sigma = 0.03; % [N/m]
    phi = 0.3;
    k = 1e-12; % [m^2]
    dssds = -1/(1-slc);
end

% calculate the normalized value of the saturation
ss = (1-sg-slc)./(1-slc);

% % load the pc functions
p = @(S)cappress(S, g, sigma, phi, k, slc, L);
dp = @(S)diffcappress(S, g, sigma, phi, k, slc, L, dssds);

% calculate the value of the relperm based on the size of the saturation
% variable
if isstruct(s)==1
    m = size(s.xvalue);
    d = ndims(s.xvalue);
    if (min(m) == 1) % 1D: only x value
        pc.xvalue = p(ss.xvalue);
        dpcds.xvalue = dp(ss.xvalue);
    elseif (d == 2) % 2D: x and y values
        pc.xvalue = p(ss.xvalue);
        pc.yvalue = p(ss.yvalue);
        dpcds.xvalue = dp(ss.xvalue);
        dpcds.yvalue = dp(ss.yvalue);
    else % 3D:x, y, and z values
        pc.xvalue = p(ss.xvalue);
        pc.yvalue = p(ss.yvalue);
        pc.zvalue = p(ss.zvalue);
        dpcds.xvalue = dp(ss.xvalue);
        dpcds.yvalue = dp(ss.yvalue);
        dpcds.zvalue = dp(ss.zvalue);
    end
else
    pc = p(ss);
    dpcds = dp(ss);
end
end


% *************************************************************************
% define the capillary pressure function
% change these two lines to define a new pc function
% *************************************************************************

function r = cappress(ss_raw, g, sigma, phi, k, slc, L)
% filter the ss values
zone1 = (ss_raw<0); zone2 = (ss_raw>1);
S = ss_raw;
S(zone1) = 0;  S(zone2) = 1;
r = (g*sigma*sqrt(phi/k)*((0.5-slc)/(1-slc))^(1/L)*(S+eps).^(-1/L));
end

function r = diffcappress(ss_raw, g, sigma, phi, k, slc, L, dssds)
% filter the ss values
zone1 = (ss_raw<0); zone2 = (ss_raw>1);
S = ss_raw;
S(zone1) = 0;  S(zone2) = 1;
r = (-1/L*dssds*g*sigma*sqrt(phi/k)*((0.5-slc)/(1-slc))^(1/L)*(S+eps).^(-1/L-1));
r(zone1) =0;    r(zone2) = 0;
end