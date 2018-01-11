function y = PQSQ_fast(x, potentialFunction, k)
%PQSQ calculates the PQSQ energy potentials u(x) in accordance with formula
%(2) in paper. 
%   y is vector of PQSQ potentials for all points.
%   x is vector of one coordinate for all data points
%   intervals is vector of intervals boundaries for one coordinate. The
%       last element MUST be Inf.
%   k is number of coordinate
%   A and B are matrices of coefficients of quadratic functions

    inds = splitPointIntervalsFast1(x,potentialFunction.intervals(k,:));
    y = potentialFunction.A(k,inds)'.*(x.^2)+potentialFunction.B(k,inds)';
end