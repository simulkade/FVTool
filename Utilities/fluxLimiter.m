function FL = fluxLimiter(varargin)
% This function returns a function handle to a flux limiter of user's
% choice.
% available flux limiters are: 'CHARM', 'HCUS', 'HQUICK', 'VanLeer',
% 'VanAlbada1', 'VanAlbada2', 'MinMod', 'SUPERBEE', 'Sweby', 'Koren',
% 'QUICK', 'MC', and 'UMIST'. Default limiter is 'SUPERBEE'. See
% <http://en.wikipedia.org/wiki/Flux_limiter>

% beta constant for Sweby and Osher 
if nargin == 2
    b = varargin{2};
else
    b=1.5;
end

% find the flux limiter function
if nargin == 1
    if strcmpi(varargin{1}, 'CHARM')
        FL = @(r)((r>0).*r.*(3*r+1)./(((r+1).^2)+eps*(r==-1)));
    elseif strcmpi(varargin{1}, 'HCUS')
        FL = @(r)(1.5*(r+abs(r))./(r+2));
    elseif strcmpi(varargin{1}, 'HQUICK')
        FL = @(r)(2.0*(r+abs(r))./((r+3)+eps*(r==-3)));
    elseif strcmpi(varargin{1}, 'ospre')
        FL = @(r)((1.5*r.*(r+1))./(r.*(r+1)+1+eps*((r.*(r+1)+1)==0)));
    elseif strcmpi(varargin{1}, 'VanLeer')
        FL = @(r)((r+abs(r))./(1+abs(r)));
    elseif strcmpi(varargin{1}, 'VanAlbada1')
        FL = @(r)((r+r.^2)./(1+r.^2));
    elseif strcmpi(varargin{1}, 'VanAlbada2')
        FL = @(r)(2*r./(1+r.^2));
    elseif strcmpi(varargin{1}, 'MinMod')
        FL = @(r)((r>0).*min(r,1));
    elseif strcmpi(varargin{1}, 'SUPERBEE')
        FL = @(r)(max(0, max(min(2*r,1), min(r,2))));
    elseif strcmpi(varargin{1}, 'Osher')
        FL = @(r)(max(0, min(r,b)));
    elseif strcmpi(varargin{1}, 'Sweby')
        FL = @(r)(max(0, max(min(b*r,1), min(r,b))));
    elseif strcmpi(varargin{1}, 'smart')
        FL = @(r)(max(0, min(4,min(0.25+0.75*r, 2*r))));
    elseif strcmpi(varargin{1}, 'Koren')
        FL = @(r)(max(0, min(2*r, min((1+2*r)/3, 2))));
    elseif strcmpi(varargin{1}, 'MUSCL')
        FL = @(r)(max(0, min(2*r, min(0.5*(1+r), 2))));
    elseif strcmpi(varargin{1}, 'QUICK')
        FL = @(r)(max(0, min(2, min(2*r, (3+r)/4))));
    elseif strcmpi(varargin{1}, 'UMIST')
        FL = @(r)(max(0, min(2, min(2*r, min((1+3*r)/4, (3+r)/4)))));
    else
        warning(['The flux limiter of your choice is not available. ' ... 
            'The SUPERBEE flux limiter is used instead.']);
        FL = @(r)(max(0, max(min(2*r,1), min(r,2))));
    end
end

if nargin == 0
    FL = @(r)(max(0, max(min(2*r,1), min(r,2))));
end