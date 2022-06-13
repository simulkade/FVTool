function dpp = fnder2(pp)
%MYFNDER returns dpp, the fist derivative of the piecewise polynomial form
% pp, which is a piecewise polynomial of one order less than pp in general
% ( same as fnder(pp,1)=fnder(pp) in Matlab's spline fitting toolbox)
% obtained from https://github.com/kurtmitman/BKM_MIT/blob/master/myfnder.m
if (pp.order == 1)
    dcoefs = zeros(size(pp.coefs));
elseif (pp.order > 1)
    dorder = pp.order-1;
    expo = repmat(fliplr(1:dorder),[pp.pieces 1]);
    dcoefs = pp.coefs(:,1:dorder).*expo;
end

dpp = mkpp(pp.breaks,dcoefs);

end

