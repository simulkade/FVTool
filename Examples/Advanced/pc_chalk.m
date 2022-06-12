% JCR library for cretaceous chalk
data{3}=[0.079 0.1
0.1 0.02
0.3 0.015
0.54 0.0
0.68 -0.4
0.747 -0.7];
% fit model to data
pc1=data{3}(:,2)*1e5; % Pa
sw1=data{3}(:,1);
sw_plot=linspace(0,1,5000);
% pc_plot=interp1(sw, pc, sw_plot, 'linear', 'extrap');
% plot(sw, pc, 'o', sw_plot, pc_plot)
% plot(diff(pc_plot)./diff(sw_plot))
% log(pc)=log(pce)-(1/labda)*log((sw-sw0)/(1-sw0-sor))
% ind0=find(pc==0, 1);
% p1=polyfit(sw(1:ind0-1), log(pc(1:ind0-1)), 1)
% plot(sw(1:ind0-1), log(pc(1:ind0-1)), 'o')

% pc=@(sw, cwi, coi, awi, aoi, swc, sor)(cwi./((sw-swc)/(1-swc)).^awi+coi./((1-sw-sor)/(1-sor)).^aoi);
% swc0=0.0789;
% sor0=0.2529;
% labda=2.4;
% f=@(x, sw)pc(sw, x(1), x(2), 1/labda, 1/labda, swc0, sor0);
% f2=@(x, sw)pc(sw, x(1), x(2), x(3), x(4), swc0, sor0);
% fw=@(x, sw)pc(sw, x(1), 0.0 , x(3), x(4), swc0, sor0);
% fo=@(x, sw)pc(sw, 0.0, x(2), x(3), x(4), swc0, sor0);
% % x=lsqcurvefit(f2, [1e3, -1e3, 1.0, 1.0],sw1, pc1)
% x=lsqnonlin(@(x)(sum((f2(x, sw1)-pc1).^2)), [0.5e3, -1e3, 1.0, 1.0])
% fit a spline
% first, do an extrapolation to assign values to sw=0 and sw=1
sw_end = [0,1];
pc_end = interp1(sw1, pc1, sw_end, 'linear', 'extrap');
pc_pp = makima([sw_end(1); sw1; sw_end(2)], [pc_end(1); pc1; pc_end(2)]);
pc = @(sw)ppval(pc_pp, sw);
pcder = @(sw)ppval(fnder(pc_pp, 1), sw);
plot(sw1, pc1, 'o', sw_plot, pc(sw_plot), '--' )%, ...
% sw_plot, fw(x, sw_plot), '--', ...
% sw_plot, fo(x, sw_plot), '-.'); 
grid;

% Note:
% Looking at the above data, it is better to go with interpolation of the
% data
% note 2: fnder requires spline toolbox
% see this: https://se.mathworks.com/matlabcentral/answers/95194-how-do-i-find-the-derivative-of-a-spline-curve-in-matlab-7-9-r2009b