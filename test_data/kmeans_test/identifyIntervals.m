function inds = identifyIntervals(x,intervals)
%identifyIntervals returns the n-by-1 vector with numbers of intervals for
%data x.
%   x is n-by-1 vector with non negative data values
%   intervals is 1-by-p vector with interval boundaries.
%Output n-by-1 vector inds contains numbers of intervals which contains
%corresponding values:
%inds(k) = q if intervals(q)<=x(k)<intervals(q+1)

[~, inds] = histc(abs(x),intervals);