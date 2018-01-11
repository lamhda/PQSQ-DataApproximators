function inds = splitPointIntervalsFast1(x,intervals)

[~, inds] = histc(abs(x),intervals);