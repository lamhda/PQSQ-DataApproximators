function y = PQSQ(x, potentialFunction, k)
%PQSQ calculates the PQSQ energy potentials u(x) in accordance with formula
%(2) in paper for kth coordinate
%   y is vector of PQSQ potentials for all points.
%   x is vector of one coordinate for all data points
%   potentialFunction is structure which define PQSQ potential function
%       (see definePotentialFunction.m for details)
%   k is number of coordinate

    inds = identifyIntervals(x,potentialFunction.intervals(k,:));
    y = potentialFunction.A(k,inds)'.*(x.^2)+potentialFunction.B(k,inds)';
end